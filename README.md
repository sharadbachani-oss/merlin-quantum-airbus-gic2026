# Merlin Quantum — GIC 2026 · Airbus — Predictive Aerodynamic Modelling

Phase-1 concept proposal, 2026 Global Quantum + AI Challenge.

| | |
|---|---|
| **Submitted proposal** | [`AIRBUS_PROPOSAL_v9.md`](AIRBUS_PROPOSAL_v9.md) — rendered as `report.pdf` |
| **Team** | Merlin Quantum, the quantum applications division of Merlin Digital (Dubai) |
| **Independent validation** | GIC 2026 **dual-track finalist** — Mitsubishi/AIST materials and QCi tracks; same framework, same instrument |
| **Repository** | https://github.com/sharadbachani-oss/merlin-quantum-airbus-gic2026 |
| **Superseded material** | `archive/` — earlier drafts and planning notes, kept for provenance; not part of the submission |

## Claim → receipt

Every measurement the proposal reports, with the file in this repository that backs it. Each entry resolves inside this repo.

| # | Measurement | Receipt |
|---:|---|---|
| 1 | This team then built the strongest classical attack on that same benchmark and measured it: one exact FFT, 7.1 s at 4,096² to 1.4×10⁻⁷ | [`floor_airbus_fft_exact.json`](floor_airbus_fft_exact.json) |
| 2 | In memory it holds at every rung unconditionally: 93–124 logical qubits against 137 GB to 35 TB of field | [`airbus_crossing_summary.json`](airbus_crossing_summary.json) |
| 3 | Benchmark parameters (statement §6), unchanged in every run | [`schrod_tgv.py`](schrod_tgv.py) |
| 4 | ii) Exact sparse Krylov solve of the transport operator (scipy.sparse.linalg.expm_multiply, measured to N = 512, law fitted on N ≥ 256) and a pseudo-spectral RK2 stepper with exact integrating factor on the identical stencil… | [`airbus_pseudospectral_comparator.json`](airbus_pseudospectral_comparator.json) |
| 5 | Appendix A): hardware waveforms correlate with the exact statevector evolution at 0.9855 and 0.920 after one damping envelope | [`qlga_result_20260810_122411.json`](qlga_result_20260810_122411.json) · [`qlga_advantage.json`](qlga_advantage.json) |
| 6 | floor_airbus_fft_exact.py (FFT solve to N = 4,096, measured), airbus_noncirculant_crossing.py (Krylov to N = 512), airbus_pseudospectral_comparator.py (time-stepper to N = 1,024); cost laws fitted on the measured rungs,… | [`floor_airbus_fft_exact.py`](floor_airbus_fft_exact.py) · [`airbus_noncirculant_crossing.py`](airbus_noncirculant_crossing.py) · [`airbus_pseudospectral_comparator.py`](airbus_pseudospectral_comparator.py) |
| 7 | schrod_tgv_results.json, constants declared; recomputed in Phase 2 with compiled-circuit constants | [`schrod_tgv_results.json`](schrod_tgv_results.json) |
| 8 | 32-core workstation; the Re ladder and the dense benchmark are receipted (airbus_bench_workstation.json) | [`airbus_bench_workstation.json`](airbus_bench_workstation.json) |
| 9 | Annex, measured, on Re-scaled grids. The tensor-train finite-volume solver on its own merits, resolution rising with Reynolds number as §4.1 requires — Re 10 at N = 64, Re 100 at N = 256, Re 1,000 at N = 1,024, Re 10⁴ at N = 4,096 | [`qtt_ladder.json`](qtt_ladder.json) · [`floor_airbus_spectral.json`](floor_airbus_spectral.json) |
| 10 | Real-time evolution of a Taylor–Green-prepared state, 42 qubits, 9 steps × 2 modes, 8,192 shots; waveform correlation 0.9855 / 0.920 | [`qlga_counts_20260810_122411.json`](raw_counts/qlga_counts_20260810_122411.json) · [`qlga_design.json`](qlga_design.json) |
| 11 | Exact FFT solve of the benchmark operator to N = 4,096 (the classical adversary this team built) | [`FLOOR_AIRBUS_FFT_EXACT.md`](FLOOR_AIRBUS_FFT_EXACT.md) |
| 12 | Krylov solve of the transport operator to N = 512 and the crossover model | [`airbus_noncirculant_crossing.json`](airbus_noncirculant_crossing.json) · [`AIRBUS_CROSSING_NONCIRCULANT.md`](AIRBUS_CROSSING_NONCIRCULANT.md) |
| 13 | Tensor-train Re ladder (N, steps, χ, bytes, L2, wall) | [`rerun_Re100_20260914.json`](rerun_Re100_20260914.json) |

## Verifying

```
python verify.py
```

Replays the headline numbers from archived counts — no credentials, no network, numpy only.
