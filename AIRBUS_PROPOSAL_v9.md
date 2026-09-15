# Resolution that stops costing what it costs today: a fault-tolerant linear-transport solver for aerodynamic design, validated on the convecting Taylor–Green vortex

Global Quantum + AI Challenge 2026 · Airbus Enterprise Challenge · Phase 1 Concept Proposal · Team Merlin Digital

Seven guideline sections · every prior-work number carries a receipt · public repository: https://github.com/sharadbachani-oss/merlin-quantum-airbus-gic2026

**Independently assessed.** This computation framework is a dual-track finalist in the GIC 2026 international quantum challenge — the MIT/Mitsubishi materials track and the QCi track — advanced by a technical panel on the operator-level method applied here.

---

## 1. Problem framing

Airbus sizes primary structure and predicts cruise fuel burn at the edge of the operating envelope, where modelling capability runs out and wind-tunnel campaigns take over. The binding constraint is degrees of freedom: on the best classical route every halving of grid spacing costs 8 to 16× in machine time and 4× in memory. On the fault-tolerant route proposed here the same refinement costs a factor of 4 in time and four qubits in memory — so design iterations affordable inside a fixed budget grow with every level added. **Where the advantage sits:** the brief asks for "a potential advantage over classical techniques by simulating a canonical fluid dynamics test case"; the convecting Taylor–Green vortex is delivered as that test case and as the calibrating instrument, the advantage one operator step away on the same solver.

**Why this benchmark admits that claim.** The 2D convecting Taylor–Green vortex is an exact Navier–Stokes solution whose nonlinear term is a pure pressure gradient, so its velocity obeys the linear advection–diffusion operator exactly (Appendix B) — the one class where fault-tolerant quantum algorithms carry a provable, not heuristic, scaling. On the periodic cell grid that operator is also circulant, giving an exact one-FFT classical solve: the referee every one of the sixteen validation rungs carries.

**The wall, measured.** The crossing: transport of a perturbation by a *non-uniform* base flow, A = −U(x,y)∂x − V(x,y)∂y + ν∇² — the linear operator of perturbation transport and linearised convection on a frozen mean flow in an aerodynamic design loop. No transform diagonalises it: the circulant structure that makes the benchmark trivially solvable is gone. Both classical routes are measured on that operator by this team, not assumed — exact sparse Krylov growing as N⁴, a pseudo-spectral stepper as N³ log N (§2) — and the better of the two is the comparator at every rung, the brief's state-of-the-art HPC comparator charged on top as a thousand-fold parallel speed-up of it.

**The instrument, calibrated — and the adversary, built.** The solver is validated at sixteen rungs against the analytic vortex, matching the exact classical solve of the same operator to 10⁻⁸–10⁻⁵ on 16–29 simulated qubits — the regime where exact classical truth exists to check against. This team then built the strongest classical attack on that same benchmark and measured it: one exact FFT, 7.1 s at 4,096² to 1.4×10⁻⁷ (`floor_airbus_fft_exact.json`). That measurement closes time-to-solution on the canonical vortex for every algorithm, and it is exactly what locates the advantage on the non-circulant transport operator, where no transform diagonalises the matrix and the attack we built has nowhere to go. Validated instrument plus self-built adversary is what licenses this solver's reading past the wall.

**The crossing, costed.** The Schrödingerised route grows as N²·polylog N on declared constants, so cores buy the classical side a fixed factor: they move the crossing, they do not close it. On one core it sits at N ≈ 69,000 per axis at Re 10⁴ on 93 logical qubits, and N ≈ 5×10⁶ at Re 100 (bounds in §8). In memory it holds at every rung unconditionally: 93–124 logical qubits against 137 GB to 35 TB of field (`airbus_crossing_summary.json`). An advantage of exponent — measured on the classical side, modelled with every quantum constant declared, and reduced in Phase 2 to one compiled number against a bar declared today (§5).

## 2. Technical approach

**Paradigm.** Fault-tolerant gate-model Schrödingerisation of the discretised transport operator, simulated exactly today at 16–29 qubits, with a quantics tensor-train finite-volume solver as the classical-hardware annex on identical stencils.

**Discretisation.** Cell averages on an N×N periodic grid, central fluxes, five-point viscous stencil: du/dt = A u with A = −U_c D_x + νL for the benchmark and A = −U D_x − V D_y + νL for the transport extension. Every comparator in the package uses this matrix.

**Benchmark parameters (statement §6), unchanged in every run** (`schrod_tgv.py`): domain length L = 2π, vortex velocity V0 = 1.0, convection velocities U_c = 1.0 and V_c = 0.0, density ρ = 1.0, background reference pressure p0 = 0.0, and viscosity ν = V0L/Re = 2π/Re — the statement's convention, **not** 1/Re, so Re 10⁴ here is ν = 6.28×10⁻⁴. ρ and p0 enter only kinematically: the solver is velocity-only, the pressure being the gradient absorbed in Appendix B.

**Schrödingerisation.** Split A = H1 + iH2 with H1 = (A + AT)/2 (dissipation, non-positive) and H2 = (A − AT)/2i (advection, Hermitian). The warped state w(p,t) = e^{−p}u(t) obeys, in the Fourier dual ξ of the auxiliary coordinate p, the unitary evolution i∂tw = (ξ⊗H1 − I⊗H2)w. Register: 2·log2N grid qubits plus n_p for p (set by the diffusive stiffness 8νT/dx² at dp = 0.05, reported per rung) plus, for the transport extension, a coefficient register of 2·log2N. Recovery from any slice p ≥ 0; readout of three design observables (kinetic-energy decay, two modal amplitudes) by amplitude estimation to ε, charged at O(1/ε).

**Cost model, constants declared.** Qubitization query count ‖H‖T + log(1/ε); five shift terms per axis, each an n-qubit incrementer of 2n² CNOT-scale gates; the coefficient oracle of the transport extension charged at twice the incrementer cost; 1 MHz logical gate rate; 10 GFLOP/s per classical core. One constant set produces both results.

**Classical comparators, all named.** (i) Exact FFT solve of the circulant benchmark operator, measured. (ii) Exact sparse Krylov solve of the transport operator (`scipy.sparse.linalg.expm_multiply`, measured to N = 512, law fitted on N ≥ 256) and a pseudo-spectral RK2 stepper with exact integrating factor on the identical stencil (measured to N = 1,024, agreeing with Krylov to 2×10⁻⁶; `airbus_pseudospectral_comparator.json`); the faster of the two binds at each rung. (iii) Dense explicit finite volume, the industrial baseline. (iv) Pseudo-spectral storage floor for the one-mode vortex, 128 bytes.

**Hardware pathway, measured.** Real-time evolution of a Taylor–Green-prepared state on 42 qubits of ibm_marrakesh (8,192 shots, nine time steps; Appendix A): hardware waveforms correlate with the exact statevector evolution at 0.9855 and 0.920 after one damping envelope (`qlga_result_20260810_122411.json`), beating a single-exponential decay hypothesis by 0.13–0.14 (`qlga_advantage.json`) — the evolution primitive the algorithm consumes, already flown.

**Annex, tensor train in finite volume.** The same stencils as quantics tensor trains: exact rank-2 shift operators, tensor-train conjugate-gradient projection, RK2; the Re-scaled ladder of §4, rank ≤ 8, L2 4×10⁻⁵ vs the analytic solution, cross-validated against dense to 10⁻⁸.

## 3. Feasibility and resource requirements

| resource | status |
|---|---|
| Solver and validation ladder | `schrod_tgv.py`, laptop, minutes: exact simulation of the algorithm's output at N = 16–128, Re 10 to 10⁴, sixteen rungs |
| Classical walls and crossing | `floor_airbus_fft_exact.py` (FFT solve to N = 4,096, measured), `airbus_noncirculant_crossing.py` (Krylov to N = 512), `airbus_pseudospectral_comparator.py` (time-stepper to N = 1,024); cost laws fitted on the measured rungs, best-of-two in `airbus_crossing_summary.json` |
| Resource model | `schrod_tgv_results.json`, constants declared; recomputed in Phase 2 with compiled-circuit constants |
| Logical register and circuit class | 93–136 logical qubits carrying exact shallow preparation, a constant-depth stencil and a phase-channel readout (below); no fault-tolerant hardware is required in Phase 2 |
| Annex compute | 32-core workstation; the Re ladder and the dense benchmark are receipted (`airbus_bench_workstation.json`) |
| Data | none; the analytic vortex and the exact classical solves are the references |

**The register this algorithm actually needs.** This is not deep variational search: the initial state is a product of sinusoids prepared exactly in O(log N) gates; the evolution is a constant-depth five-point stencil of shift operators repeated in time; the three design observables are closed-form, read through amplitude estimation — a phase/spectral quantity, not a depth-sensitive expectation sum, no optimiser loop, no ansatz search. The crossing therefore asks for 93–136 error-corrected qubits carrying one circuit family, never the full field — a register this framework reaches by its own route, not by waiting on a general fault-tolerance roadmap — and Phase 2 fixes its cost by compiling it and measuring T-count and depth. The base flow is frozen (linear class); nonlinear flows leave it through Carleman or lattice-Boltzmann embeddings (§8).

## 4. Expected impact

**The benchmark, delivered as asked (T = 0.5; L2 vs the analytic vortex; receipts in `schrod_tgv_results.json`).**

| Re | N | qubits (grid + p) | L2, Schrödingerised | L2, classical solve of the same operator | Schrödingerised vs classical |
|---|---|---|---|---|---|
| 10 | 32 | 21 | 1.41×10⁻³ | 1.43×10⁻³ | 3.5×10⁻⁵ |
| 10 | 128 | 29 | 8.03×10⁻⁵ | 8.95×10⁻⁵ | 2.1×10⁻⁵ |
| 100 | 32 | 19 | 2.13×10⁻³ | 2.13×10⁻³ | 3.5×10⁻⁶ |
| 100 | 128 | 26 | 1.33×10⁻⁴ | 1.34×10⁻⁴ | 8.2×10⁻⁶ |
| 1,000 | 128 | 23 | 1.41×10⁻⁴ | 1.41×10⁻⁴ | 4.6×10⁻⁶ |
| 10,000 | 32 | 18 | 2.27×10⁻³ | 2.27×10⁻³ | 4.4×10⁻⁸ |
| 10,000 | 128 | 22 | 1.42×10⁻⁴ | 1.42×10⁻⁴ | 1.5×10⁻⁷ |

Across all sixteen rungs (Re 10–10⁴, N 16–128) the output matches the exact classical solve to 10⁻⁸–10⁻⁵ and the analytic vortex to the discretisation error, which falls as N⁻² and is independent of Re at fixed N. **The ladder holds N fixed across Re by design:** matched-discretisation calibration; §4.1's Re-scaled resolution is carried by the annex below. **Scaling on the benchmark:** quantum memory 2·log2N + n_p qubits (45 at N = 4,096; 61–71 at 65,536) against 4N²·8 bytes of grid storage, or 128 bytes in the spectral representation of this one-mode flow; time-to-solution: the exact classical FFT solve takes 7.1 s at N = 4,096 and ≈ 39 min at N = 65,536 on one core, against 17.1 h and 319 days for the modelled fault-tolerant run at Re 10⁴ (`floor_airbus_fft_exact.json`).

**The crossing, measured on the classical side (transport by a non-uniform base flow, U = U_c − ½cos x sin y, V = ½ sin x cos y; `airbus_crossing_summary.json`).**

| Re | N | Krylov, exact | pseudo-spectral stepper | Schrödingerised at 1 MHz (oracle and readout charged) | quantum / best classical | logical qubits |
|---|---|---|---|---|---|---|
| 10,000 | 512 | 2.2 s (measured) | 25.3 s | 1.1×10³ s | 509 | 51 |
| 10,000 | 4,096 | 1.4×10³ s | 1.8×10⁴ s | 1.2×10⁵ s | 87 | 69 |
| 10,000 | 16,384 | 2.3×10⁵ s | 1.4×10⁶ s | 2.6×10⁶ s | 11 | 81 |
| 10,000 | 65,536 | 5.1×10⁷ s | 1.0×10⁸ s | 5.5×10⁷ s | 1.1 | 93 |
| 10,000 | 262,144 | 1.3×10¹⁰ s | 7.2×10⁹ s | 1.1×10⁹ s | 0.15 | 105 |
| 10,000 | 1,048,576 | 3.2×10¹² s | 5.1×10¹¹ s | 2.2×10¹⁰ s | 0.043 | 117 |
| 100 | 512 | 23.3 s (measured) | 26.2 s | 1.1×10⁵ s | 4742 | 58 |
| 100 | 4,096 | 7.4×10⁴ s | 1.8×10⁴ s | 1.2×10⁷ s | 672 | 76 |
| 100 | 16,384 | 1.9×10⁷ s | 1.4×10⁶ s | 2.7×10⁸ s | 194 | 88 |
| 100 | 65,536 | 4.8×10⁹ s | 1.0×10⁸ s | 5.5×10⁹ s | 55 | 100 |
| 100 | 262,144 | 1.2×10¹² s | 7.2×10⁹ s | 1.1×10¹¹ s | 15 | 112 |
| 100 | 1,048,576 | 3.2×10¹⁴ s | 5.1×10¹¹ s | 2.2×10¹² s | 4.3 | 124 |

![Figure 1 — Left: L2 against the analytic vortex for the Schrödingerised solver (colour) and the classical solve of the same operator (gray), Re 10–10⁴. Right: time-to-solution on the non-circulant transport operator at Re 10⁴ and 100 — measured Krylov (squares) and pseudo-spectral stepper (triangles) with their fitted cost laws (dashed) against the modelled fault-tolerant run (open circles); the dotted lines mark the single-core crossovers, ≈ 69,000 at Re 10⁴ and ≈ 5×10⁶ at Re 100.](C:/quantum ai 2026/airbus/airbus_v9_fig.png)

**Where this sits in the field.** The end-to-end quantum lattice-Boltzmann analysis (Jennings et al., arXiv:2512.03758) bounds the *nonlinear* speed-up at O(Re^{3D/8}), polynomial and conditional, with nonlinearity, readout and state preparation as the obstacles — all three removed by the linear transport class, which is why the crossing is placed there; no quantum CFD run has crossed a classical solver on time-to-solution. The quantum-inspired direction (arXiv:2512.07615, arXiv:2608.26995, PRR 7, 013112) delivers compression, not a crossing.

**Where the crossing is decisive:** perturbation transport and linearised convection on non-uniform mean flows beyond 69,000 cells per axis — where the better of our two measured classical routes needs 10⁷–10¹² s on a core and terabytes of field, the register stays under 110 logical qubits, and the design loop today simply stops.

**Annex, measured, on Re-scaled grids.** The tensor-train finite-volume solver on its own merits, resolution rising with Reynolds number as §4.1 requires — Re 10 at N = 64, Re 100 at N = 256, Re 1,000 at N = 1,024, Re 10⁴ at N = 4,096 (`airbus_bench_workstation.json`, `qtt_ladder.json`): memory ratio to the dense field 56× to 92,183×; wall-clock slower than dense at Re ≤ 1,000 and 5.2× faster at Re 10⁴ on the same grid at a 2,000× accuracy penalty; against the spectral floor of the one-mode vortex it stores 18–105× more (`floor_airbus_spectral.json`).

**Scalability to industrial relevance.** The register grows by two grid and two coefficient qubits per doubling of resolution, the base flow entering through the coefficient oracle without changing the circuit family. Three dimensions add log2N qubits per register and raise the classical exponents by one — which is where the gap widens.

**Business value, in design iterations.** The table prices the 4× against 8–16× ratio of §1: 87× at 4,096², 0.15× at 262,144², 0.043× at 1,048,576² — a 2,000-fold swing in affordable iterations over five refinement levels at Re 10⁴; in absolutes at 262,144², three design observables in 35 years of 1 MHz machine time against 229 years on one core.

## 5. Validation plan

Phase 2 is one costed three-month experiment whose output is a number. Metrics and challengers declared in advance: (1) compile the Hamiltonian-simulation circuit (Trotter and qubitization variants) with the coefficient oracle for N = 16 and report T-count and depth, replacing modelled constants with compiled ones and recomputing both crossovers; (2) extend the validation ladder to N = 256 and Re 10⁴, error and n_p per rung; (3) re-measure both classical comparators on the transport operator at N = 4,096 on the workstation, add a semi-implicit multigrid stepper as a third, and take both benchmark comparators to N = 16,384 at matched error, time and memory per rung; (4) re-fly the 42-qubit primitive with the compiled stencil circuits at N = 4 on IBM Heron: the hardware execution of one Schrödingerised step. **Pass/fail, declared now:** compiled constants hold the Re 10⁴ transport crossover inside one octave of N ≈ 69,000 — 34,500 to 138,000 per axis — at ≤ 120 logical qubits, and validation error at every rung equals the finite-volume discretisation error within a factor of two. Outside it we publish the number it lands on and the revised crossing; both outcomes are committed. The physics is established, the instrument is calibrated against exact classical truth, and the crossing point is measured by our own adversary; the remaining variable is hardware access at the scale the crossing needs — exactly what a Phase-2 PoC sprint supplies.

## 6. Hybrid / cross-domain integration

The quantum step is one stage of a classical CFD pipeline: mesh, operator assembly and base flow classical; state preparation from the analytic profile; Hamiltonian simulation on the device; amplitude estimation returns the design observables to the classical post-processor, which never reconstructs the full field. The tensor-train annex holds the same operator classically, so both routes share stencils and referees, the FFT and Krylov solves.

## 7. Team capability

Merlin Quantum is the quantum division of Merlin Digital (50+ technology FTE): Suhail Bachani (Founder & CEO, Principal Investigator), Dr. Hiro Bachani PhD (Program Director), Rohit Bachani (co-founder), Mitul Sawlani (engineering, Purdue), Mohamed Jafrun (systems), Roshan Ghazali (data), Dr. Ana Sara (validation), Gavin Vaz (platform). Prior work: 42-qubit real-time evolution on IBM Heron with gates declared before the run; the quantics tensor-train finite-volume solver and its Re ladder; 128- and 156-qubit dynamics on two devices in the sister tracks.

## 8. Scope, with treatment

(i) The canonical vortex is circulant, and our own FFT adversary is the faster route there at every Re — the measurement that placed the advantage on the non-circulant transport operator. (ii) The transport crossover is a resource-model crossing: declared quantum constants against the better of two classical cost laws this team measured, Reynolds- and hardware-dependent. Single core, N ≈ 69,000 per axis at Re 10⁴ (65,536 is the nearest tabulated power of two, ratio 1.08) and N ≈ 5×10⁶ at Re 100; against a 1,000-core node, ≈ 6×10⁷ and ≈ 8.7×10⁹. Phase 2 replaces the modelled constants with compiled ones (§5). (iii) We measured the spectral storage floor too — 128 bytes for the one-mode vortex — so the memory advantage is priced against grid storage; the transport operator's field admits no such representation. (iv) The 42-qubit job executes the evolution primitive, not a flow solve. (v) At low Re the diffusive stiffness raises the p register (n_p up to 15 at N = 128, Re 10), reported per rung. (vi) Nonlinear flows carry the polynomial bound of the brief's [9]; this crossing sits in the linear class, where the scaling is provable.

---

### Appendix A — Hardware job register and receipts

| item | device | job ID | receipt |
|---|---|---|---|
| Real-time evolution of a Taylor–Green-prepared state, 42 qubits, 9 steps × 2 modes, 8,192 shots; waveform correlation 0.9855 / 0.920 | ibm_marrakesh | `d9solsgpdb6s73e6a6q0` | `qlga_result_20260810_122411.json`, `qlga_counts_20260810_122411.json`, design `qlga_design.json` |
| Scout job | ibm_marrakesh | `d9solbntfhrs73dt1t40` | `qlga_state.json` |
| Schrödingerised solver validation ladder and resource model | CPU | — | `schrod_tgv.py`, `schrod_tgv_results.json`, `schrod_tgv_fig.png` |
| Exact FFT solve of the benchmark operator to N = 4,096 (the classical adversary this team built) | CPU | — | `floor_airbus_fft_exact.py`, `floor_airbus_fft_exact.json`, `FLOOR_AIRBUS_FFT_EXACT.md` |
| Krylov solve of the transport operator to N = 512 and the crossover model | CPU | — | `airbus_noncirculant_crossing.py`, `airbus_noncirculant_crossing.json`, `AIRBUS_CROSSING_NONCIRCULANT.md` |
| Pseudo-spectral time-stepper on the transport operator to N = 1,024 (agrees with Krylov to 2×10⁻⁶); best-of-two summary | CPU | — | `airbus_pseudospectral_comparator.py`, `airbus_pseudospectral_comparator.json`, `airbus_crossing_summary.json` |
| Tensor-train Re ladder (N, steps, χ, bytes, L2, wall) | CPU | — | `qtt_ladder.json`; Re 100 re-run `results/rerun_Re100_20260914.json` |
| Dense vs tensor-train wall-clock and error, Re 10 → 10⁴ | CPU, workstation | — | `airbus_bench_workstation.json` |
| Spectral storage floor on the canonical flow | CPU | — | `floor_airbus_spectral.json` |

### Appendix B — Why the linear operator is exact for the benchmark, and why it is circulant

For u = U_c − cos(x − U_c t) sin y·e^{−2νt}, v = sin(x − U_c t) cos y·e^{−2νt}, the convective term (u·∇)u equals −∇(¼e^{−4νt}(cos 2(x − U_c t) + cos 2y)) + U_c ∂xu: a gradient absorbed by the pressure plus mean-flow advection. What remains is u_t + U_c u_x = ν∇²u per component. On the periodic cell grid with central fluxes its matrix A = −U_c D_x + νL is a sum of shift operators with constant coefficients, hence circulant in both axes and diagonal in the two-dimensional discrete Fourier basis: exp(AT)u0 is one forward transform, one diagonal multiply and one inverse transform. The package's validation column is computed in that basis, and the floor receipt reproduces it to four digits. A non-uniform base flow multiplies the shift operators by a position-dependent field, breaks the circulant structure and removes the transform route; that is where §4's crossing sits.

### Appendix C — Claim ledger (measured · modelled · comparator · quantum attribution · cost · acceptance)

| claim | status | classical comparator | quantum attribution | cost charged | acceptance |
|---|---|---|---|---|---|
| Solver reproduces the analytic vortex to discretisation error, sixteen rungs | measured (exact simulation, 16–29 qubits) | exact classical solve of the same operator; analytic solution | algorithm output, simulated | minutes CPU | error within 2× of FV at every rung |
| Exact FFT adversary, built and measured by this team — the referee that calibrates the ladder and locates the crossing | measured | exact FFT solve, 7.1 s at N = 4,096 | — | seconds CPU | measured; frontier located |
| Memory 2·log2N + n_p qubits vs grid storage | exact | grid storage; spectral floor named | register size | — | — |
| Transport-operator crossover, fitted, at N ≈ 69,000 (Re 10⁴) and ≈ 5×10⁶ (Re 100), single core | classical side measured to N = 512 (Krylov) and 1,024 (stepper) and fitted; quantum side modelled, constants declared | best of sparse Krylov and pseudo-spectral stepper | Hamiltonian simulation with coefficient oracle | readout O(1/ε) and oracle ×2 charged | compiled constants move the Re 10⁴ crossover < 1 octave at ≤ 120 logical qubits (Phase 2) |
| 42-qubit evolution primitive on hardware | measured | exact statevector of the medium | hardware execution | 1 banked job | correlation ≥ 0.9 after one envelope (met) |
| Tensor-train annex: Re ladder, memory ratio vs dense | measured | dense FV; spectral floor | quantum-inspired | workstation hours | L2 ≤ 10⁻⁴ per rung (met) |
