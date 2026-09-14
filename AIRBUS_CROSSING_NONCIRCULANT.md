# Airbus — where the fault-tolerant route has a wall to cross (laptop, 2026-09-14 16:20)

Companion to `FLOOR_AIRBUS_FFT_EXACT.md`. Script `airbus_noncirculant_crossing.py`, receipt
`airbus_noncirculant_crossing.json`.

**Operator.** Transport of a perturbation by a non-uniform base flow on the identical periodic FV grid:
A(U) = −U(x,y)D_x − V(x,y)D_y + νL with U = U_c − ½cos x sin y, V = ½ sin x cos y (frozen TGV shape). This is the
linear operator of perturbation transport / linearised convection on a mean flow, the first extension of the
benchmark an aerodynamic design loop actually uses. No FFT diagonalises it.

**Classical comparator, measured.** Exact sparse Krylov (scipy `expm_multiply`), one laptop process:

| Re | N | Krylov wall | ‖A‖T | quantum (v8 constants, oracle ×2) at 1 MHz | logical qubits |
|---:|---:|---:|---:|---:|---:|
| 100 | 512 | 46.9 s | 1,669 | 1.1×10⁴ s | 58 |
| 10⁴ | 512 | 5.4 s | 74 | 114 s | 51 |

Cost law fitted on N ≥ 256: wall ∝ nnz·‖A‖T (∝ N⁴ at fixed Re through the diffusive norm). Extrapolated:

| Re | N | Krylov | quantum at 1 MHz | quantum / classical | logical qubits |
|---:|---:|---:|---:|---:|---:|
| 10⁴ | 4,096 | 2.7×10³ s | 1.2×10⁴ s | 4.5 | 69 |
| 10⁴ | 16,384 | 4.5×10⁵ s | 2.7×10⁵ s | 0.59 | 81 |
| 10⁴ | 65,536 | 9.8×10⁷ s | 5.5×10⁶ s | 0.056 | 93 |
| 100 | 16,384 | 3.6×10⁷ s | 2.7×10⁷ s | 0.73 | 88 |
| 100 | 65,536 | 9.3×10⁹ s | 5.5×10⁸ s | 0.060 | 100 |

**Reading.** On the non-circulant operator the crossover against a single classical core sits at N ≈ 16,384 per
axis at both Reynolds numbers, 81–88 logical qubits, and widens polynomially (N² in the model: quantum ∝ N²·polylog,
Krylov ∝ N⁴). Against a 1,000-core node the crossing moves to ≈ 2¹⁸ per axis; that must be stated. The memory side is
unconditional: 93–100 logical qubits against 8.6–137 GB.

**What this licenses in a proposal.** (1) The canonical TGV is delivered as the statement asks (validation ladder,
error/memory/time vs Re) with the one-line truth that its operator is circulant and classically exact by FFT, so no
time-to-solution advantage exists on it and none is claimed. (2) The advantage is placed on the non-circulant transport
operator: exponential in memory, polynomial (N²) in time, crossover grid and qubit count stated with declared constants,
single-core and HPC-node comparators both named. (3) Phase-2 acceptance: compiled T-counts move the crossover by less than
one octave in N; the Krylov comparator is re-measured on the workstation at N = 4,096.

## Correction, 2026-09-14 evening

Two corrections after the note above. (1) The Krylov timings were re-measured on an unloaded machine (the first run shared the laptop with an XGBoost ladder): N = 512 now 23.3 s / 2.2 s at Re 100 / 10⁴, and the Krylov-only crossover moves from 16,384 to 65,536. (2) Krylov is not the fastest classical route: a pseudo-spectral RK2 stepper with exact integrating factor (`airbus_pseudospectral_comparator.py`, measured to N = 1,024, agrees with Krylov to 2×10⁻⁶) scales as N³ log N and is the binding comparator at low Reynolds number. Best of the two (`airbus_crossing_summary.json`): single-core crossover **65,536 per axis at Re 10⁴ (93 logical qubits)**, **≈ 10⁶ at Re 100 (124 qubits)**; beyond 10⁶ against a 1,000-core node. The proposal now carries these numbers.
