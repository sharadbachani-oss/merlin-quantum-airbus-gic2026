# -*- coding: utf-8 -*-
"""airbus_noncirculant_crossing.py -- where a fault-tolerant linear-PDE solver has a wall to cross (laptop, 2026-09-14).

The canonical convecting TGV obeys A = -Uc Dx + nu L, circulant, so exp(AT)u0 is one FFT pair (floor_airbus_fft_exact).
The first Airbus-relevant extension of the same solver is transport of a perturbation by a NON-UNIFORM base flow,
A(U) = -U(x,y) Dx - V(x,y) Dy + nu L, the linear operator of every perturbation / passive-scalar / linearised-convection
problem on a frozen mean flow.  No FFT diagonalises it.  The classical exact route is sparse Krylov (scipy expm_multiply),
measured here; the quantum route is the same Schrodingerisation with a coefficient oracle, modelled with the v8 constants
plus a declared oracle factor.  Base flow: U = Uc - a cos x sin y, V = a sin x cos y (frozen TGV shape, a = 0.5).
"""
import json, math, time, os, sys
import numpy as np
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import expm_multiply
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import schrod_tgv as S

A_BASE = 0.5
def op_noncirc(N, nu):
    dx = S.LBOX / N; xs = (np.arange(N) + 0.5) * dx
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    U = S.UC - A_BASE * np.cos(X) * np.sin(Y); V = A_BASE * np.sin(X) * np.cos(Y)
    e = np.ones(N); Sh = diags([e[:-1], e[:1]], [1, -(N - 1)], shape=(N, N))
    D1 = (Sh - Sh.T) / (2 * dx); L1 = (Sh + Sh.T - 2 * identity(N)) / dx ** 2; I = identity(N)
    Dx = kron(D1, I, format="csr"); Dy = kron(I, D1, format="csr"); L = kron(L1, I) + kron(I, L1)
    A = (-diags(U.ravel()) @ Dx - diags(V.ravel()) @ Dy + nu * L).tocsr()
    return A, X, Y

def classical_krylov(N, nu, T):
    A, X, Y = op_noncirc(N, nu)
    u0 = (-np.cos(X) * np.sin(Y)).ravel()          # perturbation initial field (TGV shape, mean removed)
    t0 = time.perf_counter(); uT = expm_multiply(A, u0, start=0.0, stop=T, num=2, endpoint=True)[-1]
    wall = time.perf_counter() - t0
    normA = float(abs(A).sum(axis=1).max())         # inf-norm bound
    return wall, normA, A.nnz, float(np.linalg.norm(uT) / np.linalg.norm(u0))

def quantum_model(N, Re, T, oracle_factor=2.0):
    """v8 resources() with the coefficient oracle charged: each stencil term multiplied by U(x,y) needs a block-encoding
    of a smooth function of 2n qubits; declared as oracle_factor x the incrementer cost (arithmetic of the same depth)."""
    r = S.resources(N, Re, T)
    return dict(logical_qubits=r["logical_qubits"] + 2 * int(math.log2(N)),   # + coefficient register
                quantum_gates=r["quantum_gates"] * oracle_factor, quantum_s_at_1MHz=r["quantum_gates"] * oracle_factor / 1e6)

if __name__ == "__main__":
    T = 0.5; rows = []
    for Re in (100, 10000):
        nu = S.V0 * S.LBOX / Re
        for N in (64, 128, 256, 512):
            w, nA, nnz, decay = classical_krylov(N, nu, T)
            q = quantum_model(N, Re, T)
            rows.append(dict(Re=Re, N=N, krylov_wall_s=w, normA=nA, nnz=nnz, normAT=nA * T, decay=decay, **q))
            print(f"Re {Re:6d} N {N:4d}: Krylov {w:8.2f} s  |A|T {nA*T:9.1f}  nnz {nnz:9d} | quantum {q['quantum_gates']:.2e} gates "
                  f"= {q['quantum_s_at_1MHz']:.3g} s @1MHz, {q['logical_qubits']} logical qubits", flush=True)
    # classical cost law: wall ~ c * nnz * |A|T  (Taylor/Krylov steps x sparse matvec); fit c on N>=256, extrapolate
    fit = [(r["krylov_wall_s"] / (r["nnz"] * r["normAT"])) for r in rows if r["N"] >= 256]
    c = float(np.median(fit))
    ext = []
    for Re in (100, 10000):
        nu = S.V0 * S.LBOX / Re
        for N in (1024, 4096, 16384, 65536):
            dx = S.LBOX / N; nA = 8 * nu / dx ** 2 + (S.UC + A_BASE) / dx * 2; nnz = 5 * N * N
            w = c * nnz * nA * T; q = quantum_model(N, Re, T)
            ext.append(dict(Re=Re, N=N, krylov_wall_s_extrapolated=w, quantum_s_at_1MHz=q["quantum_s_at_1MHz"],
                            logical_qubits=q["logical_qubits"], quantum_over_krylov=q["quantum_s_at_1MHz"] / w,
                            classical_bytes=16 * N * N * 8))
            print(f"extrap Re {Re:6d} N {N:6d}: Krylov {w:.3g} s | quantum {q['quantum_s_at_1MHz']:.3g} s | q/c {q['quantum_s_at_1MHz']/w:.2e} | {q['logical_qubits']} qubits")
    xo = {}
    for Re in (100, 10000):
        cand = [e for e in ext if e["Re"] == Re and e["quantum_over_krylov"] < 1]
        xo[str(Re)] = min(cand, key=lambda e: e["N"])["N"] if cand else None
    out = dict(card="non-circulant transport operator: measured Krylov vs modelled Schrödingerisation",
               base_flow=f"U = Uc - {A_BASE} cos x sin y, V = {A_BASE} sin x cos y (frozen TGV shape)", T=T,
               classical_fit_seconds_per_nnz_per_normAT=c, measured=rows, extrapolated=ext, crossover_N_per_Re=xo,
               assumptions=dict(v8_constants="schrod_tgv.resources (1 MHz logical, eps 1e-4, 3 observables, amplitude estimation)",
                                oracle_factor=2.0, classical_rate="this laptop, single process scipy expm_multiply",
                                note="classical comparator is the exact Krylov solve, not explicit stepping; no FFT route exists for variable coefficients"))
    json.dump(out, open("airbus_noncirculant_crossing.json", "w"), indent=1)
    print("crossover N per Re:", xo, "-> airbus_noncirculant_crossing.json")
