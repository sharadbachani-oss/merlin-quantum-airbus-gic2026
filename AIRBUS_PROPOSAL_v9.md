# A fault-tolerant linear-transport solver, validated on the convecting Taylor–Green vortex, with its classical wall measured and its crossing placed

Global Quantum + AI Challenge 2026 · Airbus Enterprise Challenge · Phase 1 Concept Proposal · Team Merlin Digital (GIC 2026 dual-track finalist)

Seven guideline sections · every prior-work number carries a receipt · public repository: https://github.com/sharadbachani-oss/merlin-quantum-airbus-gic2026

---

## 1. Problem framing

Aerodynamic design is bounded by the cost of solving the flow equations at the resolutions extreme operating conditions demand, and the statement asks for a quantum solver graded on the canonical case: the 2D convecting Taylor–Green vortex, exact solution known, with time-to-solution, memory and error quantified as the Reynolds number rises, and with the Reynolds number taken "to the maximum required to demonstrate possible advantage".

Two facts about that benchmark decide where a credible advantage can and cannot sit, and this proposal is built on both. **First,** the convecting Taylor–Green vortex is an exact solution of the full nonlinear Navier–Stokes equations whose nonlinear term is a pure pressure gradient, so its velocity obeys the linear advection–diffusion operator exactly (Appendix B). Linear evolution is the one class where fault-tolerant quantum algorithms carry a provable, not heuristic, scaling: the N×N field is held in 2·log₂N qubits and evolved by Hamiltonian simulation. **Second,** on the periodic cell grid that operator is circulant, so the exact classical solution is one Fourier transform pair with no time stepping: we measured it at 6.5 s for a 4,096² grid to 1.4×10⁻⁷ error on a laptop (`floor_airbus_fft_exact.json`). No quantum algorithm beats that on the canonical vortex, and we do not claim one; a proposal that claimed a time-to-solution crossover on this benchmark against explicit time stepping would be measuring against the wrong classical algorithm, and we say so because our own earlier draft did exactly that before the floor caught it.

**Credible advantage, stated and bounded.** The crossing sits one step from the benchmark, on the same solver: transport of a perturbation by a *non-uniform* base flow, A = −U(x,y)∂ₓ − V(x,y)∂ᵧ + ν∇², the linear operator of every perturbation, passive-scalar and linearised-convection problem on a frozen mean flow in an aerodynamic design loop. No Fourier transform diagonalises it. Two classical routes are measured on it, not assumed: exact sparse Krylov exponentiation, whose cost grows as N⁴ at fixed Reynolds number, and a pseudo-spectral explicit time-stepper, whose cost grows as N³ log N; the better of the two is the comparator at every rung. The Schrödingerised route grows as N²·polylog N with the same declared constants. Against the better classical route on a single core the crossover is at N ≈ 69,000 per axis at Re 10⁴ (93 logical qubits) and near 5×10⁶ per axis at Re 100 (136 logical qubits), where the diffusive stiffness inflates the quantum cost; against a thousand-core node it lies near 6×107 per axis at Re 10⁴ and beyond any grid in use at Re 100, and we say so. In memory the register is 93–124 logical qubits against 137 GB to 35 TB, unconditionally (`airbus_crossing_summary.json`). That is an advantage of exponent, measured on the classical side and modelled with every constant on the table on the quantum side, and it is what a successful PoC makes concrete.

## 2. Technical approach

**Paradigm.** Fault-tolerant gate-model algorithm (Schrödingerisation of the discretised transport operator), simulated exactly today at 16–29 qubits; a quantics tensor-train finite-volume solver as the classical-hardware annex on identical stencils; the real-time evolution primitive flown at 42 qubits on IBM Heron as the hardware pathway.

**Discretisation.** Cell averages on an N×N periodic grid, central fluxes, five-point viscous stencil: du/dt = A u with A = −U_c D_x + νL for the benchmark and A = −U D_x − V D_y + νL for the transport extension. The dense, tensor-train, FFT and Krylov comparators in the package use the same matrix, so every comparison is on one operator against one analytic or exact solution.

**Schrödingerisation.** Split A = H₁ + iH₂ with H₁ = (A + Aᵀ)/2 (dissipation, non-positive) and H₂ = (A − Aᵀ)/2i (advection, Hermitian). The warped state w(p,t) = e^{−p}u(t) obeys, in the Fourier dual ξ of the auxiliary coordinate p, the unitary evolution i∂ₜŵ = (ξ⊗H₁ − I⊗H₂)ŵ. Register: 2·log₂N grid qubits plus n_p for p (set by the diffusive stiffness 8νT/dx² at dp = 0.05, reported per rung) plus, for the transport extension, a coefficient register of 2·log₂N. Recovery from any slice p ≥ 0; readout of three design observables (kinetic-energy decay, two modal amplitudes) by amplitude estimation to ε, charged at O(1/ε).

**Cost model, constants declared.** Qubitization query count ‖H‖T + log(1/ε); five shift terms per axis, each an n-qubit incrementer of 2n² CNOT-scale gates; the coefficient oracle of the transport extension charged at twice the incrementer cost; 1 MHz logical gate rate; 10 GFLOP/s per classical core. Nothing is fitted to the outcome; the same constants produce the benchmark's "no crossover" and the transport operator's crossover.

**Classical comparators, all named.** (i) Exact FFT solve of the circulant benchmark operator (the correct comparator on the canonical vortex, measured). (ii) Exact sparse Krylov solve of the transport operator (`scipy.sparse.linalg.expm_multiply`, measured to N = 512, cost law fitted on N ≥ 256) and a pseudo-spectral RK2 time-stepper with exact integrating factor for diffusion on the identical stencil (measured to N = 1,024, agreeing with Krylov to 2×10⁻⁶; `airbus_pseudospectral_comparator.json`); the faster of the two is the comparator at each rung. (iii) Dense explicit finite volume on identical stencils (in package; named as the industrial baseline, not as the comparator). (iv) Pseudo-spectral storage floor for the one-mode vortex, 128 bytes.

**Hardware pathway, measured.** Real-time evolution of a Taylor–Green-prepared state on 42 qubits of ibm_marrakesh (one banked job, 8,192 shots, nine time steps, two modes): hardware waveforms correlate with the exact statevector evolution at 0.993 and 0.971 after one damping envelope, and a detuned hypothesis loses by 0.14. It is the time-stepping primitive the algorithm consumes, run on today's device; it is not a flow solve and is not claimed as one.

**Annex, tensor train in finite volume.** The same stencils as quantics tensor trains: exact rank-2 shift operators, tensor-train conjugate-gradient projection, RK2; Re 10 → 10⁴, N to 4,096, rank ≤ 8, L2 4×10⁻⁵ vs the analytic solution, cross-validated against dense to 10⁻⁸; graded on its own merits in §4.

## 3. Feasibility and resource requirements

| resource | status |
|---|---|
| Solver and validation ladder | `schrod_tgv.py`, laptop, minutes: exact simulation of the algorithm's output at N = 16–128, Re 10 to 10⁴, sixteen rungs |
| Classical walls and crossing | `floor_airbus_fft_exact.py` (FFT solve to N = 4,096, measured), `airbus_noncirculant_crossing.py` (Krylov to N = 512), `airbus_pseudospectral_comparator.py` (time-stepper to N = 1,024); cost laws fitted on the measured rungs, best-of-two in `airbus_crossing_summary.json` |
| Resource model | `schrod_tgv_results.json`, constants declared; recomputed in Phase 2 with compiled-circuit constants |
| Fault-tolerant hardware | not required in Phase 2; the crossing sits at 93–136 logical qubits, inside published fault-tolerant roadmaps; the primitive is already flown at 42 physical qubits |
| Annex compute | 32-core workstation; the Re ladder and the dense benchmark are receipted (`airbus_bench_workstation.json`) |
| Data | none; the analytic vortex and the exact classical solves are the references |

Assumptions: the Taylor–Green initial state is a product of sinusoids prepared in O(log N) gates; readout targets three observables, never the full field; the transport base flow is frozen (linear class). Constraint: nonlinear flows leave the linear class and enter the polynomial bound of the brief's [9] through Carleman or lattice-Boltzmann embeddings; Phase 2 states where that bound bites and does not claim past it.

## 4. Expected impact

**The benchmark, delivered as asked (T = 0.5; L2 vs the analytic vortex; receipts in `schrod_tgv_results.json`).**

| Re | N | qubits (grid + p) | L2, Schrödingerised | L2, classical solve of the same operator | Schrödingerised vs classical |
|---|---|---|---|---|---|
| 100 | 32 | 19 | 2.13×10⁻³ | 2.13×10⁻³ | 3.5×10⁻⁶ |
| 100 | 128 | 26 | 1.33×10⁻⁴ | 1.34×10⁻⁴ | 8.2×10⁻⁶ |
| 1,000 | 128 | 23 | 1.41×10⁻⁴ | 1.41×10⁻⁴ | 4.6×10⁻⁶ |
| 10,000 | 32 | 18 | 2.27×10⁻³ | 2.27×10⁻³ | 4.4×10⁻⁸ |
| 10,000 | 128 | 22 | 1.42×10⁻⁴ | 1.42×10⁻⁴ | 1.5×10⁻⁷ |

Across all sixteen rungs (Re 10–10⁴, N 16–128) the algorithm's output matches the exact classical solve of the same operator to 10⁻⁸–10⁻⁵ and the analytic vortex to the discretisation error, which falls as N⁻² and is independent of Re at fixed N. **Scaling on the benchmark:** error N⁻²; quantum memory 2·log₂N + n_p qubits (45 at N = 4,096; 61–71 at 65,536) against 4N²·8 bytes of grid storage, or 128 bytes in the spectral representation of this one-mode flow; time-to-solution: the exact classical FFT solve takes 6.5 s at N = 4,096 and ≈ 37 min at N = 65,536 on one laptop core, against 1.7 h and 32 days for the modelled fault-tolerant run at Re 10⁴. **No time-to-solution advantage exists on the canonical vortex, at any Reynolds number, and none is claimed.**

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

**The field's direction, and its wall.** The end-to-end quantum lattice-Boltzmann analysis (Jennings et al., arXiv:2512.03758) bounds the speed-up for nonlinear flows at O(Re^{3D/8}), polynomial and conditional, with nonlinearity, readout and state preparation as the obstacles; no quantum CFD run has crossed a classical solver on time-to-solution. The quantum-inspired direction — tensor-network lattice Boltzmann (arXiv:2512.07615), tensor-train compressible turbulence (arXiv:2608.26995), GPU quantics turbulence (PRR 7, 013112) — delivers compression, not a crossing. On the canonical vortex the wall is absolute: the operator is circulant and the exact classical solve is one FFT.

**The crossing.** We take the benchmark at its word and the field at its wall. The solver is built and validated on the analytic vortex at sixteen rungs; the classical comparator that actually beats it on the benchmark is named and measured; and the advantage is placed where the same solver meets an operator no transform diagonalises, with the crossover grid, the qubit count and both classical comparators stated as numbers. **Where it is decisive:** perturbation transport and linearised convection on non-uniform mean flows at high Reynolds number and grids beyond 69,000 per axis, where the best classical route needs 10⁷–10¹² s on a core and terabytes of field, and the register stays under 110 logical qubits.

**Annex, measured.** The tensor-train finite-volume solver on its own merits: memory ratio to the dense field 56× to 92,183× across Re 10 → 10⁴; wall-clock slower than dense at Re ≤ 1,000 and 5.2× faster at Re 10⁴ on the same grid at a 2,000× accuracy penalty (`airbus_bench_workstation.json`); against the spectral floor of the one-mode vortex it stores 18–105× more (`floor_airbus_spectral.json`). It is the classical-hardware route to the same operator, not the advantage claim.

**Scalability to industrial relevance.** The register grows by two grid qubits and two coefficient qubits per doubling of resolution; the operator is a five-point stencil at any N; the base flow enters through the coefficient oracle with no change to the circuit family; three dimensions add log₂N qubits per register and raise the classical exponents by one, which is where the exponent gap widens. The nonlinear extension is Carleman linearisation into the same Schrödingerisation, inside the polynomial bound of the brief's [9].

**Business value, bounded.** In Airbus's units the value is design iterations per wind-tunnel campaign. As a scenario computed from the table, not a consequence of the model: a linearised-transport case at 262,144² on a non-uniform mean flow at Re 10⁴ returns its three design observables in 3.5 years on a 1 MHz fault-tolerant machine against 23 years on one classical core or ≈ 8 days on a 1,000-core node; the fault-tolerant route is the one whose cost rises by 4× per halving of grid spacing while the best classical one rises by 8–16×, and that ratio, not today's constant, is what decides the fleet-scale case.

## 5. Validation plan

Phase 2, three months, metrics and challengers declared in advance: (1) compile the Hamiltonian-simulation circuit (Trotter and qubitization variants) with the coefficient oracle for N = 16 and report T-count and depth, replacing modelled constants with compiled ones and recomputing both crossovers; (2) extend the validation ladder to N = 256 and Re 10⁴, error and n_p per rung; (3) re-measure both classical comparators on the transport operator at N = 4,096 on the workstation, add a semi-implicit multigrid stepper as a third, and the FFT comparator on the benchmark to N = 16,384, so every classical cost law carries two more octaves of measurement; (4) pseudo-spectral comparator on the benchmark at matched error, time and memory per rung; (5) re-fly the 42-qubit primitive with the compiled stencil circuits at N = 4 on IBM Heron as the hardware execution of one Schrödingerised step. **Acceptance:** compiled constants move the Re 10⁴ transport crossover by less than one octave in N with ≤ 120 logical qubits; validation error at every rung equals the finite-volume discretisation error within a factor of two; either outcome is reported.

## 6. Hybrid / cross-domain integration

The quantum step is one stage of a classical CFD pipeline: mesh, operator assembly and the base-flow field classical; initial-state preparation from the analytic profile; Hamiltonian simulation on the fault-tolerant device; amplitude estimation returns the design observables to the classical post-processor, which never reconstructs the full field. The tensor-train annex is the same operator held classically, so the two routes share stencils, verification and referee; the FFT and Krylov solves are the referees for both.

## 7. Team capability

Merlin Quantum is the quantum division of Merlin Digital (50+ technology FTE): Suhail Bachani (Founder & CEO, Principal Investigator), Dr. Hiro Bachani PhD (Program Director), Rohit Bachani (co-founder), Mitul Sawlani (engineering, Purdue), Mohamed Jafrun (systems), Roshan Ghazali (data), Dr. Ana Sara (validation), Gavin Vaz (platform). Prior work: 42-qubit real-time evolution on IBM Heron with gates declared before the run; the quantics tensor-train finite-volume solver and its Re ladder; 128- and 156-qubit real-time dynamics on two devices in the sister tracks; GIC 2026 dual-track finalist, with the same discipline of measuring the classical wall before placing a claim.

## 8. Scope, with treatment

(i) On the canonical vortex there is no time-to-solution advantage against the exact FFT solve, which we measured; the benchmark is delivered for validation and scaling, and the advantage is placed on the non-circulant transport operator. (ii) The transport crossover is a resource-model crossing with declared constants against the better of two measured classical cost laws; it is Reynolds-dependent (65,536 per axis at Re 10⁴, ≈ 10⁶ at Re 100) and against a 1,000-core node it lies beyond 10⁶ per axis; Phase 2 replaces the constants with compiled ones and adds a third classical route. (iii) The memory advantage on the one-mode vortex is against grid storage; the spectral representation holds that flow in 128 bytes, and the transport operator's field has no such representation. (iv) The 42-qubit job executes the evolution primitive, not a flow solve. (v) At low Re the diffusive stiffness raises the p register (n_p up to 15 at N = 128, Re 10), reported per rung. (vi) Nonlinear flows enter the polynomial bound of the brief's [9]; nothing here claims past it. (vii) The tensor-train annex is faster than dense only at matched grid at Re 10⁴, not at matched error.

---

### Appendix A — Hardware job register and receipts

| item | device | job ID | receipt |
|---|---|---|---|
| Real-time evolution of a Taylor–Green-prepared state, 42 qubits, 9 steps × 2 modes, 8,192 shots; waveform correlation 0.993 / 0.971 | ibm_marrakesh | `d9solsgpdb6s73e6a6q0` | `qlga_result_20260810_122411.json`, `qlga_counts_20260810_122411.json`, design `qlga_design.json` |
| Scout job | ibm_marrakesh | `d9solbntfhrs73dt1t40` | `qlga_state.json` |
| Schrödingerised solver validation ladder and resource model | CPU | — | `schrod_tgv.py`, `schrod_tgv_results.json`, `schrod_tgv_fig.png` |
| Exact FFT solve of the benchmark operator to N = 4,096 (the classical wall on the canonical vortex) | CPU | — | `floor_airbus_fft_exact.py`, `floor_airbus_fft_exact.json`, `FLOOR_AIRBUS_FFT_EXACT.md` |
| Krylov solve of the transport operator to N = 512 and the crossover model | CPU | — | `airbus_noncirculant_crossing.py`, `airbus_noncirculant_crossing.json`, `AIRBUS_CROSSING_NONCIRCULANT.md` |
| Pseudo-spectral time-stepper on the transport operator to N = 1,024 (agrees with Krylov to 2×10⁻⁶); best-of-two summary | CPU | — | `airbus_pseudospectral_comparator.py`, `airbus_pseudospectral_comparator.json`, `airbus_crossing_summary.json` |
| Tensor-train Re ladder (N, steps, χ, bytes, L2, wall) | CPU | — | `qtt_ladder.json`; Re 100 re-run `results/rerun_Re100_20260914.json` |
| Dense vs tensor-train wall-clock and error, Re 10 → 10⁴ | CPU, workstation | — | `airbus_bench_workstation.json` |
| Spectral storage floor on the canonical flow | CPU | — | `floor_airbus_spectral.json` |

### Appendix B — Why the linear operator is exact for the benchmark, and why it is circulant

For u = U_c − cos(x − U_c t) sin y·e^{−2νt}, v = sin(x − U_c t) cos y·e^{−2νt}, the convective term (u·∇)u equals −∇(¼e^{−4νt}(cos 2(x − U_c t) + cos 2y)) + U_c ∂ₓu: a gradient absorbed by the pressure plus mean-flow advection. What remains is u_t + U_c u_x = ν∇²u per component. On the periodic cell grid with central fluxes its matrix A = −U_c D_x + νL is a sum of shift operators with constant coefficients, hence circulant in both axes and diagonal in the two-dimensional discrete Fourier basis: exp(AT)u₀ is one forward transform, one diagonal multiply and one inverse transform. The Schrödingerised solver in the package simulates its own output in exactly that basis, which is how the package's validation column is computed; the floor receipt reproduces that column to four digits. A non-uniform base flow multiplies the shift operators by a position-dependent field, breaks the circulant structure, and removes the transform route; that is where §4's crossing sits.

### Appendix C — Claim ledger (measured · modelled · comparator · quantum attribution · cost · acceptance)

| claim | status | classical comparator | quantum attribution | cost charged | acceptance |
|---|---|---|---|---|---|
| Solver reproduces the analytic vortex to discretisation error, sixteen rungs | measured (exact simulation, 16–29 qubits) | exact classical solve of the same operator; analytic solution | algorithm output, simulated | minutes CPU | error within 2× of FV at every rung |
| No time-to-solution advantage on the canonical vortex | measured | exact FFT solve, 6.5 s at N = 4,096 | — | seconds CPU | stated, not contested |
| Memory 2·log₂N + n_p qubits vs grid storage | exact | grid storage; spectral floor named | register size | — | — |
| Transport-operator crossover at N ≈ 65,536 (Re 10⁴) and ≈ 10⁶ (Re 100), single core | classical side measured to N = 512 (Krylov) and 1,024 (stepper) and fitted; quantum side modelled, constants declared | best of sparse Krylov and pseudo-spectral stepper | Hamiltonian simulation with coefficient oracle | readout O(1/ε) and oracle ×2 charged | compiled constants move the Re 10⁴ crossover < 1 octave at ≤ 120 logical qubits (Phase 2) |
| 42-qubit evolution primitive on hardware | measured | exact statevector of the medium | hardware execution | 1 banked job | correlation ≥ 0.9 after one envelope (met) |
| Tensor-train annex: Re ladder, memory ratio vs dense | measured | dense FV; spectral floor | quantum-inspired | workstation hours | L2 ≤ 10⁻⁴ per rung (met) |
