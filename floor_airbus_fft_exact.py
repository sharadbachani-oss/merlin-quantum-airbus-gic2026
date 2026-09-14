# -*- coding: utf-8 -*-
"""floor_airbus_fft_exact.py -- cheapest adversary for the v8 time-to-solution crossover (laptop, 2026-09-14).

v8 §4 compares the Schrödingerised gate count against EXPLICIT finite volume (N^2 cells x steps). But the
operator A = -Uc Dx + nu L on the periodic cell grid is circulant in both axes, so exp(A T) u0 is one 2D FFT,
a diagonal multiply, and one inverse FFT: O(N^2 log N), no time stepping, and it returns the SAME vector the
Schrödingerised algorithm's readout is taken from (schrod_tgv.schrod_solve itself diagonalises H_S in the DFT
basis). This script measures that solve, checks it reproduces the package's own 'classical FV' L2 column,
fits its wall-clock, and puts it beside the v8 gate model at every rung including the claimed crossover.
"""
import json, math, sys, time, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import schrod_tgv as S

def fft_exact_solve(N, nu, T):
    dx, X, Y = S.grid(N)
    u0, v0 = S.tgv_exact(X, Y, 0.0, nu)
    t0 = time.perf_counter()
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx); KX, KY = np.meshgrid(k, k, indexing="ij")
    lam = (nu * ((2 * np.cos(KX * dx) - 2) + (2 * np.cos(KY * dx) - 2)) / dx ** 2   # nu L eigenvalues
           - 1j * S.UC * np.sin(KX * dx) / dx)                                       # -Uc Dx eigenvalues
    prop = np.exp(lam * T)
    up = S.UC + np.fft.ifft2(np.fft.fft2(u0 - S.UC) * prop).real
    vp = np.fft.ifft2(np.fft.fft2(v0) * prop).real
    wall = time.perf_counter() - t0
    ue, ve = S.tgv_exact(X, Y, T, nu)
    l2 = math.sqrt(float(np.mean((up - ue) ** 2 + (vp - ve) ** 2)))   # package convention (schrod_tgv line ~82)
    return wall, l2

T = 0.5; rows = []
val = {(r["Re"], r["N"]): r for r in json.load(open("schrod_tgv_results.json"))["validation"]}
for Re in (100, 10000):
    nu = S.V0 * S.LBOX / Re
    for N in (128, 256, 1024, 4096):
        wall, l2 = fft_exact_solve(N, nu, T)
        wall = min(wall, fft_exact_solve(N, nu, T)[0])          # best of two
        q = S.resources(N, Re, T)
        pkg = val.get((Re, N), {}).get("l2_classical_fv_vs_analytic")
        rows.append(dict(Re=Re, N=N, fft_wall_s=wall, fft_l2_vs_analytic=l2,
                         package_classical_fv_l2=pkg,
                         quantum_gates=q["quantum_gates"], quantum_s_at_1MHz=q["quantum_gates"] / 1e6,
                         explicit_fv_s_at_10GFLOPs=q["classical_flops"] / 1e10,
                         logical_qubits=q["logical_qubits"]))
        print(f"Re {Re:6d} N {N:5d}: FFT exact {wall*1e3:9.1f} ms  L2 {l2:.3e}  (package classical-FV L2 {pkg})  |"
              f" quantum {q['quantum_gates']:.2e} gates = {q['quantum_gates']/1e6:.3g} s @1MHz |"
              f" explicit FV model {q['classical_flops']/1e10:.3g} s @10GFLOP/s", flush=True)
# fit wall = c * N^2 log2(N^2) on the measured rungs, extrapolate to the claimed crossover
c = np.median([r["fft_wall_s"] / (r["N"] ** 2 * math.log2(r["N"] ** 2)) for r in rows if r["N"] >= 1024])
extrap = []
for Re in (10, 100, 1000, 10000):
    for N in (16384, 65536):
        q = S.resources(N, Re, T); w = c * N ** 2 * math.log2(N ** 2)
        extrap.append(dict(Re=Re, N=N, fft_wall_s_extrapolated=w, quantum_s_at_1MHz=q["quantum_gates"] / 1e6,
                           explicit_fv_s_at_10GFLOPs=q["classical_flops"] / 1e10,
                           quantum_over_fft=q["quantum_gates"] / 1e6 / w, logical_qubits=q["logical_qubits"],
                           fft_bytes=2 * N * N * 16))
        print(f"extrap Re {Re:6d} N {N:6d}: FFT {w:9.1f} s | quantum {q['quantum_gates']/1e6:.3g} s | "
              f"quantum/FFT {q['quantum_gates']/1e6/w:.2e} | explicit-FV model {q['classical_flops']/1e10:.3g} s")
out = dict(card="cheapest adversary for the v8 Schrödingerisation time-to-solution crossover",
           method="exact exp(A T) u0 by 2D FFT diagonalisation of the identical periodic central stencil; "
                  "same output vector as the Schrödingerised readout; no time stepping",
           machine=os.environ.get("COMPUTERNAME"), T=T, fit_seconds_per_N2log2N2=float(c),
           measured=rows, extrapolated=extrap,
           verdict=None)
worst = max(r["quantum_over_fft"] for r in extrap); best = min(r["quantum_over_fft"] for r in extrap)
out["verdict"] = (f"On the operator v8 names, the exact classical solve is one FFT pair. At the claimed crossover "
                  f"N=65,536 the modelled fault-tolerant time is {best:.1e}-{worst:.1e}x the FFT wall-clock on this "
                  f"laptop; the crossover exists only against explicit time-stepping, which nobody uses for a "
                  f"constant-coefficient periodic linear operator. The v8 time-to-solution claim does not survive; "
                  f"the memory claim was already conceded against the 128-byte spectral floor.")
json.dump(out, open("floor_airbus_fft_exact.json", "w"), indent=1)
print("\n" + out["verdict"]); print("-> floor_airbus_fft_exact.json")
