"""Schrodingerization solver for the 2D convecting Taylor-Green vortex (FTQC track), simulated exactly.

The convecting TGV  u = Uc - cos(x-Uc t) sin(y) e^{-2 nu t},  v = sin(x-Uc t) cos(y) e^{-2 nu t}  (Vc = 0)
is an exact Navier-Stokes solution whose nonlinear term is a pure pressure gradient, so the velocity obeys the LINEAR
advection-diffusion operator exactly:   u_t + Uc u_x = nu Lap u  (same for v).  On an N x N periodic finite-volume grid
with central fluxes this is  du/dt = A u,  A = -Uc Dx + nu L  (real, non-normal, dissipative).

Schrodingerization (Jin, Liu, Yu 2022): split A = H1 + i H2 with H1 = (A + A^T)/2 (= nu L, symmetric <= 0) and
H2 = (A - A^T)/(2i) (advection, Hermitian). Warp with an auxiliary coordinate p:  w(p, t) = e^{-p} u(t) for p > 0.
Then  w_t = -H1 w_p + i H2 w,  and in the Fourier dual xi of p the evolution is UNITARY:
       i d/dt  w_hat(xi) = ( xi H1 - H2 ) w_hat(xi)      ->      H_S = xi (x) H1  -  I (x) H2      (Hermitian)
u(t) is recovered from any slice p >= 0:  u(t) = e^{p} w(p, t).  The quantum register is (x, y, p): 2 log2 N + n_p qubits.

This script simulates the algorithm's output exactly (expm_multiply of H_S on the joint statevector), validates it
against the analytic TGV at Re = 10, 100, 1000 and N = 16, 32, 64, separates Schrodingerization error from
discretisation error (compare to expm(A t) u0), and writes a gate-level resource model with the classical crossover.
Writes schrod_tgv_results.json and figs: schrod_tgv_error.png, schrod_tgv_resources.png.
"""
import os, sys, time, json, math
import numpy as np
from scipy.sparse import diags, identity, kron, csr_matrix
from scipy.sparse.linalg import expm_multiply, expm
HERE = os.path.dirname(os.path.abspath(__file__))
UC, V0, LBOX = 1.0, 1.0, 2 * math.pi

def grid(N):
    dx = LBOX / N; xs = (np.arange(N) + 0.5) * dx
    X, Y = np.meshgrid(xs, xs, indexing="ij"); return dx, X, Y

def tgv_exact(X, Y, t, nu):
    e = math.exp(-2 * nu * t)
    return UC - np.cos(X - UC * t) * np.sin(Y) * e, np.sin(X - UC * t) * np.cos(Y) * e

def operators(N, nu):
    """A = -Uc Dx + nu L on the N x N periodic cell grid, central differences (identical stencils to qtt_fv / dense)."""
    dx = LBOX / N
    e = np.ones(N)
    S = diags([e[:-1], e[:1]], [1, -(N - 1)], shape=(N, N))          # periodic shift  f_{i+1}
    Dx1 = (S - S.T) / (2 * dx)                                        # central first derivative
    L1 = (S + S.T - 2 * identity(N)) / dx ** 2                        # central second derivative
    I = identity(N)
    Dx = kron(Dx1, I, format="csr"); L = kron(L1, I, format="csr") + kron(I, L1, format="csr")
    A = (-UC * Dx + nu * L).tocsr()
    H1 = (nu * L).tocsr()                     # symmetric part (dissipation)
    H2 = ((A - A.T) / 2j).tocsr()             # Hermitian part of the advection (A - A^T = -2 Uc Dx, Dx antisymmetric)
    return A, H1, H2

def schrod_solve(N, nu, T):
    """Exact simulation of the discretised Schrodingerised algorithm's output state.
    On the periodic grid H1 = nu L and H2 (advection) are both circulant, hence commute and are diagonal in the 2D DFT
    basis; H_S = xi (x) H1 - I (x) H2 is therefore diagonal in (kx, ky, xi). The p register is kept DISCRETE (M points,
    dp = 0.05) so its truncation / periodisation error is measured, not assumed. Recovery from the slice p in [0.5, 3]."""
    dx, X, Y = grid(N); u0, v0 = tgv_exact(X, Y, 0.0, nu)
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx); KX, KY = np.meshgrid(k, k, indexing="ij")
    h1 = nu * ((2 * np.cos(KX * dx) - 2) + (2 * np.cos(KY * dx) - 2)) / dx ** 2      # eigenvalues of nu L   (<= 0)
    h2 = -UC * np.sin(KX * dx) / dx                                                    # eigenvalues of (A - A^T)/2i
    lam_max = 8 * nu / dx ** 2
    DP = 0.05; p_lo = -(lam_max * T + 4.0); p_hi = 4.0
    M = int(2 ** math.ceil(math.log2((p_hi - p_lo) / DP))); M = min(M, 2 ** 15); n_p = int(math.log2(M))
    ps = np.linspace(p_lo, p_hi, M, endpoint=False); dp = ps[1] - ps[0]
    xis = 2 * np.pi * np.fft.fftfreq(M, d=dp)
    wp_hat = np.fft.fft(np.exp(-np.abs(ps))) / M                                       # warp profile in the dual variable
    sel = np.where((ps >= 0.5) & (ps <= 3.0))[0]
    E_win = np.exp(1j * np.outer(ps[sel] - p_lo, xis)) * np.exp(ps[sel])[:, None]      # inverse DFT at the window (index convention: p - p_lo), times e^{p}
    out = {}; t0 = time.time()
    for name, f0 in (("u", u0 - UC), ("v", v0)):
        F0 = np.fft.fft2(f0).ravel()                                                   # grid modes
        h1f = h1.ravel(); adv = np.exp(1j * h2.ravel() * T) * F0
        w_win = np.zeros((len(sel), N * N), dtype=complex)
        CH = max(1, int(2 ** 26 // (N * N)))                                           # chunk the dual variable: <= 1 GB per block
        for c0 in range(0, M, CH):
            c1 = min(M, c0 + CH)
            w_hat_T = (wp_hat[c0:c1, None] * np.exp(-1j * np.outer(xis[c0:c1], h1f) * T)) * adv[None, :]
            w_win += E_win[:, c0:c1] @ w_hat_T                                         # (n_win, N^2) : e^{p} w(p, T)
        FT = w_win.mean(axis=0)
        rec = np.fft.ifft2(FT.reshape(N, N)).real.ravel()
        classical = np.fft.ifft2((np.exp((h1 + 1j * h2) * T).ravel() * F0).reshape(N, N)).real.ravel()   # exact classical solve of the same operator
        out[name] = dict(rec=rec, classical=classical)
    wall = time.time() - t0
    uT, vT = tgv_exact(X, Y, T, nu)
    err_schrod = math.sqrt(((out["u"]["rec"] + UC - uT.ravel()) ** 2 + (out["v"]["rec"] - vT.ravel()) ** 2).mean())
    err_classical = math.sqrt(((out["u"]["classical"] + UC - uT.ravel()) ** 2 + (out["v"]["classical"] - vT.ravel()) ** 2).mean())
    err_vs_classical = math.sqrt(((out["u"]["rec"] - out["u"]["classical"]) ** 2 + (out["v"]["rec"] - out["v"]["classical"]) ** 2).mean())
    return dict(N=N, Re=round(V0 * LBOX / nu), T=T, n_p=n_p, p_points=M, qubits=2 * int(math.log2(N)) + n_p,
                l2_vs_analytic=err_schrod, l2_classical_fv_vs_analytic=err_classical, l2_schrod_vs_classical=err_vs_classical, wall_s=wall)

# ----------------------------- resource model (fault-tolerant) -----------------------------
def resources(N, Re, T, eps=1e-4, eps_read=1e-2):
    """Gate-level estimate for simulating H_S = xi (x) H1 - I (x) H2 to time T within eps (qubitization / LCU).
    H1, H2 are sums of 2D shift operators (5-point stencils): s = 5 terms, each shift an n-qubit incrementer.
    Queries ~ ||H_S|| T + log(1/eps); each query ~ s * (incrementer cost) gates.  Classical: explicit FV N^2 cells x steps."""
    nu = V0 * LBOX / Re; dx = LBOX / N; n = int(math.log2(N))
    lam1 = 8 * nu / dx ** 2                    # ||H1||
    lam2 = UC / dx                             # ||H2||
    DP = 0.05; M = int(2 ** math.ceil(math.log2((lam1 * T + 8.0) / DP))); n_p = int(math.log2(M))   # same rule as the solver
    xi_max = math.pi / DP                      # Nyquist of the p grid
    normH = xi_max * lam1 + lam2
    queries = normH * T + math.log(1 / eps)
    incrementer = 2 * n ** 2                   # Toffoli-free-ancilla incrementer on n qubits (CNOT-count scale)
    gates_per_query = 5 * 2 * incrementer + 4 * (2 * n + n_p)   # 5 stencil terms x 2 axes + bookkeeping
    q_gates = queries * gates_per_query
    # readout: m observables (kinetic-energy decay, two modal amplitudes) to precision eps via amplitude estimation
    m_obs = 3; q_total = q_gates * m_obs / eps_read                 # readout: 3 design observables to precision eps_read (1%) by amplitude estimation, O(1/eps_read) oracle calls each, charged in full
    logical_qubits = 2 * n + n_p + 6
    # classical explicit finite volume, same stencils, same accuracy class
    dt = min(0.25 * dx / UC, 0.2 * dx ** 2 / nu); steps = T / dt
    c_flops = 2 * N ** 2 * 30 * steps          # ~30 flops per cell per step, u and v
    c_bytes = 4 * N ** 2 * 8
    return dict(N=N, Re=Re, logical_qubits=logical_qubits, quantum_gates=q_total, classical_flops=c_flops,
                classical_bytes=c_bytes, quantum_steps_equiv=queries, classical_steps=steps)

def crossover(Re, T, gate_rate=1e6, flop_rate=1e10):
    """Smallest N (power of two) where fault-tolerant time-to-solution < classical explicit FV time-to-solution."""
    for n in range(4, 40):
        N = 2 ** n; r = resources(N, Re, T)
        if r["quantum_gates"] / gate_rate < r["classical_flops"] / flop_rate: return N, r
    return None, None

if __name__ == "__main__":
    T = 0.5
    rows = []
    for Re in (10, 100, 1000, 10000):
        nu = V0 * LBOX / Re
        for N in (16, 32, 64, 128):
            r = schrod_solve(N, nu, T); rows.append(r)
            print(f"Re {Re:5d} N {N:3d} qubits {r['qubits']:2d}  L2 schrod {r['l2_vs_analytic']:.2e}  classical FV {r['l2_classical_fv_vs_analytic']:.2e}  |schrod-classical| {r['l2_schrod_vs_classical']:.2e}  {r['wall_s']:.1f}s", flush=True)
    res = []
    for Re in (10, 100, 1000, 10000):
        for N in (64, 256, 1024, 4096, 16384, 65536):
            res.append(resources(N, Re, T))
    xo = {Re: crossover(Re, T) for Re in (10, 100, 1000, 10000)}
    for Re, (N, r) in xo.items():
        print(f"crossover Re {Re}: N = {N}  (qubits {r['logical_qubits']}, quantum gates {r['quantum_gates']:.2e}, classical flops {r['classical_flops']:.2e}, classical bytes {r['classical_bytes']:.2e})" if N else f"crossover Re {Re}: none below 2^40")
    json.dump(dict(validation=rows, resources=res, crossover={str(k): (v[0], v[1]) for k, v in xo.items()},
                   assumptions=dict(T=T, eps=1e-4, gate_rate_logical_per_s=1e6, classical_flop_rate=1e10, p_register='dp=0.05 over [-(lam_max T+4), 4]',
                                    stencil="central 5-point, identical to the dense/QTT solvers", readout="3 observables via amplitude estimation")),
              open(os.path.join(HERE, "schrod_tgv_results.json"), "w"), indent=1)
    try:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
        for Re in (10, 100, 1000):
            rr = [r for r in rows if r["Re"] == Re]
            ax[0].loglog([r["N"] for r in rr], [r["l2_vs_analytic"] for r in rr], "o-", label=f"Schrödingerised, Re {Re}")
            ax[0].loglog([r["N"] for r in rr], [r["l2_classical_fv_vs_analytic"] for r in rr], "x--", color="gray", alpha=.6)
        ax[0].set_xlabel("N (grid per axis)"); ax[0].set_ylabel("velocity L2 vs analytic TGV"); ax[0].legend(fontsize=7); ax[0].set_title("error scaling (gray: classical FV, same stencil)", fontsize=8)
        for Re in (100, 10000):
            rr = [r for r in res if r["Re"] == Re]
            ax[1].loglog([r["N"] for r in rr], [r["quantum_gates"] / 1e6 for r in rr], "o-", label=f"quantum, Re {Re} (s at 1 MHz logical)")
            ax[1].loglog([r["N"] for r in rr], [r["classical_flops"] / 1e10 for r in rr], "x--", label=f"classical FV, Re {Re} (s at 10 GFLOP/s)")
        ax[1].set_xlabel("N"); ax[1].set_ylabel("time-to-solution (s)"); ax[1].legend(fontsize=7); ax[1].set_title("time-to-solution model, T = 0.5", fontsize=8)
        plt.tight_layout(); plt.savefig(os.path.join(HERE, "schrod_tgv_fig.png"), dpi=160); print("-> schrod_tgv_fig.png")
    except Exception as e: print("fig skipped:", e)
    print("-> schrod_tgv_results.json")
