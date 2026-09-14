# -*- coding: utf-8 -*-
"""airbus_pseudospectral_comparator.py -- the fastest standard classical route for the non-circulant transport operator
(laptop, 2026-09-14 evening): pseudo-spectral RK2 with exact integrating factor for diffusion; variable-coefficient
advection evaluated in real space, derivatives in Fourier space (2/3 dealiasing not needed for a linear operator).
Cost per step ~ 6 FFTs of N^2; steps ~ T / (0.25 dx / Umax) ∝ N  ->  wall ∝ N^3 log N.  Measured at N = 256..1024,
checked against the exact Krylov solve at N = 256, then extrapolated and put beside the v8/v9 quantum model.
"""
import json, math, time, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import schrod_tgv as S
from airbus_noncirculant_crossing import op_noncirc, quantum_model, A_BASE
from scipy.sparse.linalg import expm_multiply

def ps_solve(N, nu, T):
    dx = S.LBOX / N; xs = (np.arange(N) + 0.5) * dx; X, Y = np.meshgrid(xs, xs, indexing="ij")
    U = S.UC - A_BASE * np.cos(X) * np.sin(Y); V = A_BASE * np.sin(X) * np.cos(Y)
    u0 = -np.cos(X) * np.sin(Y)
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx); KX, KY = np.meshgrid(k, k, indexing="ij")
    # identical discrete stencil eigenvalues as the FV operator (central differences), so the same operator is solved
    ikx = 1j * np.sin(KX * dx) / dx; iky = 1j * np.sin(KY * dx) / dx
    lam = nu * ((2 * np.cos(KX * dx) - 2) + (2 * np.cos(KY * dx) - 2)) / dx ** 2
    umax = float(np.abs(U).max() + np.abs(V).max()); dt = 0.25 * dx / umax; steps = int(math.ceil(T / dt)); dt = T / steps
    E1 = np.exp(lam * dt); Eh = np.exp(lam * dt / 2)
    def adv(uh):
        ux = np.fft.ifft2(ikx * uh).real; uy = np.fft.ifft2(iky * uh).real
        return np.fft.fft2(-(U * ux + V * uy))
    t0 = time.perf_counter(); uh = np.fft.fft2(u0)
    for _ in range(steps):                      # RK2 (Heun) with integrating factor
        k1 = adv(uh); u1 = E1 * (uh + dt * k1); k2 = adv(u1)
        uh = E1 * uh + 0.5 * dt * (E1 * k1 + k2)
    wall = time.perf_counter() - t0
    return wall, steps, np.fft.ifft2(uh).real

T = 0.5; rows = []
for Re in (100, 10000):
    nu = S.V0 * S.LBOX / Re
    for N in (256, 512, 1024):
        w, steps, u = ps_solve(N, nu, T); w = min(w, ps_solve(N, nu, T)[0]) if N <= 512 else w
        err = None
        if N == 256:
            A, X, Y = op_noncirc(N, nu); u0 = (-np.cos(X) * np.sin(Y)).ravel()
            uk = expm_multiply(A, u0, start=0.0, stop=T, num=2, endpoint=True)[-1]
            err = float(np.sqrt(np.mean((u.ravel() - uk) ** 2)) / np.sqrt(np.mean(uk ** 2)))
        q = quantum_model(N, Re, T)
        rows.append(dict(Re=Re, N=N, ps_wall_s=w, steps=steps, rel_err_vs_krylov=err, quantum_s_at_1MHz=q["quantum_s_at_1MHz"], logical_qubits=q["logical_qubits"]))
        print(f"Re {Re:6d} N {N:5d}: pseudo-spectral {w:8.2f} s ({steps} steps){'  rel err vs Krylov %.1e' % err if err else ''} | quantum {q['quantum_s_at_1MHz']:.3g} s", flush=True)
c = float(np.median([r["ps_wall_s"] / (r["N"] ** 3 * math.log2(r["N"] ** 2)) for r in rows if r["N"] >= 512]))
ext = []
for Re in (100, 10000):
    for N in (4096, 16384, 65536, 262144, 1048576):
        q = quantum_model(N, Re, T); w = c * N ** 3 * math.log2(N ** 2)
        ext.append(dict(Re=Re, N=N, ps_wall_s_extrapolated=w, quantum_s_at_1MHz=q["quantum_s_at_1MHz"], quantum_over_ps=q["quantum_s_at_1MHz"] / w, logical_qubits=q["logical_qubits"]))
        print(f"extrap Re {Re:6d} N {N:7d}: pseudo-spectral {w:.3g} s | quantum {q['quantum_s_at_1MHz']:.3g} s | q/c {q['quantum_s_at_1MHz']/w:.2e} | {q['logical_qubits']} qubits")
xo = {str(Re): min([e["N"] for e in ext if e["Re"] == Re and e["quantum_over_ps"] < 1] or [None]) for Re in (100, 10000)}
json.dump(dict(card="pseudo-spectral explicit time-stepper vs modelled Schrödingerisation on the non-circulant transport operator",
               method="RK2 + exact integrating factor for diffusion, derivatives in Fourier space, same discrete stencil as the FV operator; dt = 0.25 dx/Umax",
               fit_seconds_per_N3log2N2=c, measured=rows, extrapolated=ext, crossover_N_per_Re_single_core=xo, T=T),
          open("airbus_pseudospectral_comparator.json", "w"), indent=1)
print("single-core crossover N per Re:", xo)
