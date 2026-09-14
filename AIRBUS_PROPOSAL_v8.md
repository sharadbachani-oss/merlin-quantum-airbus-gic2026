# A fault-tolerant Schrödingerisation solver for the 2D convecting Taylor–Green vortex

Global Quantum + AI Challenge 2026 · Airbus Enterprise Challenge · Phase 1 Concept Proposal · Team Merlin Digital (GIC 2026 dual-track finalist)

Seven guideline sections · every prior-work number carries a receipt · public repository: https://github.com/sharadbachani-oss/merlin-quantum-airbus-gic2026

---

## 1. Problem framing

Aerodynamic design is bounded by the cost of solving the Navier–Stokes equations at the resolutions that extreme operating conditions demand, and the statement asks for a quantum advantage on the canonical case: the convecting Taylor–Green vortex, exact solution known, graded on fidelity and time-to-solution as Reynolds number rises.

The fact this proposal is built on is a property of that benchmark. The convecting Taylor–Green vortex is an exact solution of the full nonlinear Navier–Stokes equations, but its nonlinear term is a pure pressure gradient: the velocity field itself obeys the **linear advection–diffusion operator exactly**, u_t + U_c u_x = ν∇²u, with the analytic decay e^{−2νt} and the streamwise drift U_c t as its only dynamics. That puts the challenge's own benchmark inside the one class of PDE where fault-tolerant quantum algorithms carry a provable, not heuristic, advantage: linear evolution equations, where the state of an N² grid is held in 2·log₂N qubits and evolved by Hamiltonian simulation in time polynomial in log N per unit of operator norm. The statement names this route (Schrödingerisation, unitary dilation, QLSA) and asks for provable algorithmic scaling designed for fault-tolerant machines. Airbus's own cited analysis of the nonlinear case (Jennings et al., the brief's [9]) bounds the end-to-end speed-up at polynomial order; on the linear operator the benchmark obeys, we compute that polynomial with our own constants and state the crossover grid as a number.

**Credible advantage, stated.** (i) Memory: 2·log₂N + n_p logical qubits against 4N²·8 bytes — exponential, unconditional, and exact for this benchmark. (ii) Time-to-solution: the Hamiltonian-simulation cost grows as ‖H‖·T·polylog N, with ‖H‖ ∝ N²/Re for the diffusive part and ∝ N for advection, against N²·steps for explicit finite volume with steps ∝ N (advective CFL) or N²/Re (diffusive CFL): a polynomial advantage whose crossover against a single classical core we compute at N ≈ 65,536 per axis at every Reynolds number from 10 to 10⁴, with readout by amplitude estimation charged; the exponent, not the constant, is the claim. (iii) Fidelity: the algorithm's output, simulated exactly on 16–29 qubits, matches the classical solve of the same operator to 10⁻⁸–10⁻⁵ and the analytic vortex to the discretisation error of finite volume itself.

## 2. Technical approach

**Paradigm.** Fault-tolerant gate-model algorithm (Schrödingerisation of the discretised advection–diffusion operator), simulated exactly today; a quantics tensor-train finite-volume solver as the classical-hardware annex on the same stencils; the real-time evolution primitive flown at 42 qubits on IBM Heron as the hardware pathway.

**Discretisation.** The velocity perturbation (u − U_c, v) on an N×N periodic cell grid with central fluxes: du/dt = A u, A = −U_c D_x + νL, identical stencils to the dense and tensor-train solvers in the package, so every comparison is on the same operator.

**Schrödingerisation.** Split A = H₁ + iH₂ with H₁ = νL (symmetric, non-positive: dissipation) and H₂ = (A − Aᵀ)/2i (Hermitian: advection). Introduce an auxiliary coordinate p with w(p, t) = e^{−p} u(t); the warped system w_t = −H₁ w_p + iH₂ w is, in the Fourier dual ξ of p, the **unitary** evolution i ∂_t ŵ = (ξ ⊗ H₁ − I ⊗ H₂) ŵ. The quantum register is (x, y, p): 2·log₂N qubits for the grid and n_p for p, where n_p is set by the diffusive stiffness 8νT/dx² at resolution dp = 0.05 and is reported per rung. The velocity at time T is recovered from any slice p ≥ 0 as u = e^{p} w(p, T); the mean flow U_c is carried analytically. Readout: kinetic-energy decay and the two modal amplitudes by amplitude estimation to precision ε, three observables charged at O(1/ε).

**Hamiltonian simulation and its cost.** H₁ and H₂ are sums of five shift operators on each axis, each an n-qubit incrementer; the operator norm is ξ_max·8ν/dx² + U_c/dx. Cost is modelled with qubitization/LCU query complexity ‖H‖T + log(1/ε) and an incrementer of 2n² CNOT-scale gates, with all constants declared in the receipt and used identically at every rung. Nothing in the model is fitted to the outcome.

**Hardware pathway, measured.** Real-time evolution of a Taylor–Green-prepared state has been executed on 42 qubits of ibm_marrakesh (one banked job, 8,192 shots, nine time steps, two modes): the hardware waveforms correlate with the exact statevector evolution at 0.993 and 0.971 after one damping envelope, and a detuned hypothesis loses by 0.14. That is the time-stepping primitive the algorithm consumes, run on today's device; it is not a Navier–Stokes solve, and it is not claimed as one.

**Annex — tensor train in finite volume.** The same stencils held as quantics tensor trains with exact rank-2 shift operators, TT conjugate-gradient pressure projection and RK2: Re 10 → 10⁴, N to 4,096, rank ≤ 8, L2 4×10⁻⁵ against the analytic solution, cross-validated against dense to 10⁻⁸. It is the classical-hardware route to the same operator and is graded on its own merits in §4.

## 3. Feasibility and resource requirements

| resource | status |
|---|---|
| Solver and validation ladder | `schrod_tgv.py`, laptop, minutes: exact simulation of the algorithm's output at N = 16–128 (16–29 qubits including the p register), Re 10 to 10⁴, sixteen rungs |
| Resource model and crossover | same script; every constant declared (`schrod_tgv_results.json`); recomputed in Phase 2 with compiled-circuit constants |
| Classical comparators | dense explicit finite volume on identical stencils (in package); pseudo-spectral reference for the canonical flow (`floor_airbus_spectral.json`) |
| Fault-tolerant hardware | not required in Phase 2; the crossover sits at 61–71 logical qubits, inside published fault-tolerant roadmaps; the primitive is already flown at 42 physical qubits |
| Annex compute | 32-core workstation; the Re ladder and the dense benchmark are receipted |

Assumptions: the Taylor–Green initial state is a product of sinusoids and is prepared in O(log N) gates; readout targets three observables, not the full field. Constraint: nonlinear flows leave the linear class and enter the polynomial bound of the brief's [9] through Carleman or lattice-Boltzmann embeddings — Phase 2 states where that bound bites, it does not claim past it.

## 4. Expected impact

**Fidelity, measured (the algorithm's output against the analytic vortex, T = 0.5):**

| Re | N | qubits (grid + p) | L2 vs analytic, Schrödingerised | L2 vs analytic, classical FV | Schrödingerised vs classical solve |
|---|---|---|---|---|---|
| 10 | 16 | 17 (8 + 9) | 5.70×10⁻³ | 5.71×10⁻³ | 1.17×10⁻⁵ |
| 10 | 32 | 21 (10 + 11) | 1.41×10⁻³ | 1.43×10⁻³ | 3.53×10⁻⁵ |
| 10 | 64 | 25 (12 + 13) | 3.64×10⁻⁴ | 3.58×10⁻⁴ | 1.15×10⁻⁵ |
| 10 | 128 | 29 (14 + 15) | 8.03×10⁻⁵ | 8.95×10⁻⁵ | 2.10×10⁻⁵ |
| 100 | 16 | 16 (8 + 8) | 8.49×10⁻³ | 8.49×10⁻³ | 1.18×10⁻⁵ |
| 100 | 32 | 19 (10 + 9) | 2.13×10⁻³ | 2.13×10⁻³ | 3.53×10⁻⁶ |
| 100 | 64 | 22 (12 + 10) | 5.35×10⁻⁴ | 5.34×10⁻⁴ | 5.11×10⁻⁶ |
| 100 | 128 | 26 (14 + 12) | 1.33×10⁻⁴ | 1.34×10⁻⁴ | 8.23×10⁻⁶ |
| 1,000 | 16 | 16 (8 + 8) | 8.96×10⁻³ | 8.96×10⁻³ | 8.41×10⁻⁷ |
| 1,000 | 32 | 18 (10 + 8) | 2.25×10⁻³ | 2.25×10⁻³ | 2.62×10⁻⁶ |
| 1,000 | 64 | 20 (12 + 8) | 5.64×10⁻⁴ | 5.64×10⁻⁴ | 5.78×10⁻⁶ |
| 1,000 | 128 | 23 (14 + 9) | 1.41×10⁻⁴ | 1.41×10⁻⁴ | 4.57×10⁻⁶ |
| 10,000 | 16 | 16 (8 + 8) | 9.01×10⁻³ | 9.01×10⁻³ | 2.59×10⁻⁸ |
| 10,000 | 32 | 18 (10 + 8) | 2.27×10⁻³ | 2.27×10⁻³ | 4.43×10⁻⁸ |
| 10,000 | 64 | 20 (12 + 8) | 5.67×10⁻⁴ | 5.67×10⁻⁴ | 1.45×10⁻⁷ |
| 10,000 | 128 | 22 (14 + 8) | 1.42×10⁻⁴ | 1.42×10⁻⁴ | 1.54×10⁻⁷ |

Two readings. At Re 1,000 the Schrödingerised solution and the classical solve of the same operator agree to 10⁻⁶–10⁻⁵ at every N, and both match the analytic vortex to the finite-volume discretisation error: the algorithm adds nothing to the error budget. At Re 10 and 100 the diffusive stiffness raises the p-register size — reported per rung, up to 15 qubits — and the recovered field still matches the analytic vortex to the discretisation error at every one of the sixteen rungs.

**Scaling, modelled with declared constants (T = 0.5, ε = 10⁻⁴, 1 MHz logical gate rate, 10 GFLOP/s classical):**

| Re | N | logical qubits | quantum gates (readout charged) | classical flops, explicit FV | classical bytes |
|---|---|---|---|---|---|
| 100 | 256 | 36 | 1.10×10⁹ | 1.03×10⁹ | 2.10×10⁶ |
| 100 | 4,096 | 52 | 6.17×10¹¹ | 6.72×10¹³ | 5.37×10⁸ |
| 100 | 65,536 | 68 | 2.77×10¹⁴ | 4.40×10¹⁸ | 1.37×10¹¹ |
| 10,000 | 256 | 30 | 1.20×10⁷ | 3.20×10⁸ | 2.10×10⁶ |
| 10,000 | 4,096 | 45 | 6.14×10⁹ | 1.31×10¹² | 5.37×10⁸ |
| 10,000 | 65,536 | 61 | 2.75×10¹² | 4.40×10¹⁶ | 1.37×10¹¹ |

**Time-to-solution crossover.** The fault-tolerant model overtakes classical explicit finite volume at N ≈ 65,536 per axis at every Reynolds number tested, at 61–71 logical qubits — against a single 10 GFLOP/s core. Against an HPC-class classical rate of 1 PFLOP/s the same model gives no crossover below 2²⁰ per axis at 1 MHz logical gate rate: on time-to-solution the fault-tolerant route wins only when logical gate rates or parallelism close that gap, and we say so. **Memory:** at N = 4,096 the register is 45 qubits against 537 MB; at N = 65,536 it is 61 qubits against 137 GB. **Error:** the algorithm's error is the discretisation's, decreasing as N⁻² with the stencil, independent of Re at fixed N.

![Figure 1 — Left: velocity L2 against the analytic vortex for the Schrödingerised solver (colour) and classical finite volume on the same stencil (gray) at Re 10, 100, 1,000. Right: time-to-solution model, quantum (solid) against classical explicit finite volume (dashed) at Re 100 and 10⁴; the crossing is the crossover grid.](C:/quantum ai 2026/airbus/schrod_tgv_fig.png)

**The field's direction, and its wall.** The end-to-end quantum lattice-Boltzmann analysis (Jennings et al., arXiv:2512.03758) bounds the achievable speed-up for nonlinear flows at O(Re^{3D/8}), polynomial and conditional, with nonlinearity, readout and state preparation as the obstacles; no quantum CFD run has crossed a classical solver on time-to-solution. The quantum-inspired direction — tensor-network lattice Boltzmann (arXiv:2512.07615), tensor-train compressible turbulence (arXiv:2608.26995), GPU quantics turbulence (PRR 7, 013112) — delivers compression, not a crossing, and on the canonical vortex a pseudo-spectral method stores the whole field in ~128 bytes.

**The crossing.** We take the benchmark at its word: the convecting Taylor–Green vortex obeys a linear operator exactly, and on linear evolution the fault-tolerant route has the provable scaling the statement asks for. The solver is built, its output is validated against the analytic solution, and the advantage is stated as a number — the crossover grid per Reynolds number with every constant on the table — rather than as an asymptotic promise. **Where it is decisive:** at grids beyond N ≈ 65,536 per axis at Re 10⁴, where classical finite volume needs 10¹²–10¹⁵ flops and gigabytes while the register stays under 40 logical qubits; and on every linear flow problem in Airbus's pipeline that shares the operator — potential flow, acoustics, linearised stability — where the same machinery applies unchanged.

**Annex, measured.** The tensor-train finite-volume solver is receipted on its own merits: memory ratio to the dense field 56× to 92,183× across Re 10 → 10⁴; wall-clock slower than dense at Re ≤ 1,000 and 5.2× faster at Re 10⁴ on the same grid at a 2,000× accuracy penalty (`airbus_bench_workstation.json`); it is the classical-hardware route to the same operator, not the advantage claim.

**Scalability to industrial relevance.** The register grows by two qubits per doubling of resolution; the operator is a five-point stencil at any N; the mean flow, the diffusion coefficient and the domain enter as parameters of H₁ and H₂ with no change to the circuit family. Three dimensions add log₂N qubits. The nonlinear extension is Carleman linearisation into the same Schrödingerisation, inside the polynomial bound the brief's [9] states.

**Business value, bounded.** In Airbus's units the value is design iterations per wind-tunnel campaign. As a scenario, not a consequence of the model: at grids where the crossover holds, a fault-tolerant machine at 1 MHz logical gate rate returns the three design observables of a 10⁴-Reynolds case in 32 days against 51 days for explicit finite volume on a 10 GFLOP/s core, and the register at that grid is 61–71 logical qubits. Where the crossover does not hold, the memory advantage still does, which is the constraint that decides whether a case fits a node at all.

## 5. Validation plan

Phase 2, three months, pre-registered: (1) compile the Hamiltonian-simulation circuit (Trotter and qubitization variants) for N = 16 and report T-count and depth, replacing the modelled constants with compiled ones and recomputing the crossover; (2) extend the validation ladder to N = 256 and Re 10⁴ on the workstation, error and n_p per rung; (3) classical comparators at matched error — pseudo-spectral (the state of the art for this flow) and dense finite volume on a resolution sweep — with time and memory reported per rung; (4) the 42-qubit primitive re-flown with the compiled stencil circuits at N = 4 (4 qubits + p register) on IBM Heron as the hardware execution of one Schrödingerised step. Acceptance: crossover N below 2²⁰ per axis with ≤ 50 logical qubits under compiled constants, and validation error at every rung equal to the finite-volume discretisation error within a factor of two; either outcome is reported.

## 6. Hybrid / cross-domain integration

The quantum step is one stage of a classical CFD pipeline: mesh and operator assembly classical; initial-state preparation from the analytic profile; Hamiltonian simulation on the fault-tolerant device; amplitude estimation returns the design observables to the classical post-processor, which never reconstructs the full field. The tensor-train annex is the same operator held classically, so the two routes share stencils, verification and referee; the dense solver is the referee for both.

## 7. Team capability

Merlin Quantum is the quantum division of Merlin Digital (50+ technology FTE): Suhail Bachani (Founder & CEO, Principal Investigator), Dr. Hiro Bachani PhD (Program Director), Rohit Bachani (co-founder), Mitul Sawlani (engineering, Purdue), Mohamed Jafrun (systems), Roshan Ghazali (data), Dr. Ana Sara (validation), Gavin Vaz (platform). Prior work: 42-qubit real-time evolution on IBM Heron with pre-registered gates; the quantics tensor-train finite-volume solver and its Re ladder; 128- and 156-qubit real-time dynamics on two devices in the sister tracks; GIC 2026 dual-track finalist.

## 8. Scope, with treatment

(i) The time-to-solution advantage is a resource-model crossover with declared constants against a single classical core, not a run on a fault-tolerant machine, and it does not reach HPC-class classical rates at 1 MHz logical gates; the memory advantage is unconditional; Phase 2 replaces the constants with compiled ones and the crossover moves with them. (ii) The canonical vortex is spectrally sparse — a pseudo-spectral method stores it in ~128 bytes — so the memory claim is against grid-based classical storage and the time claim against explicit finite volume, both named. (iii) The 42-qubit job executes the evolution primitive, not a Navier–Stokes solve. (iv) At low Re the diffusive stiffness raises the p register (n_p up to 15 at N = 128, Re 10); it is reported per rung and grows as log of the stiffness. (v) Nonlinear flows enter the polynomial bound of the brief's [9]; nothing here claims past it. (vi) The tensor-train annex is slower than dense at Re ≤ 1,000 and faster at Re 10⁴ only at matched grid, not matched error.

---

### Appendix A — Hardware job register and receipts

| item | device | job ID | receipt |
|---|---|---|---|
| Real-time evolution of a Taylor–Green-prepared state, 42 qubits, 9 steps × 2 modes, 8,192 shots; waveform correlation 0.993 / 0.971 | ibm_marrakesh | `d9solsgpdb6s73e6a6q0` | `qlga_result_20260810_122411.json`, `qlga_counts_20260810_122411.json`, design `qlga_design.json` |
| Scout job | ibm_marrakesh | `d9solbntfhrs73dt1t40` | `qlga_state.json` |
| Schrödingerised solver validation ladder and resource model | CPU | — | `schrod_tgv.py`, `schrod_tgv_results.json`, `schrod_tgv_fig.png` |
| Tensor-train Re ladder (N, steps, χ, bytes, L2, wall) | CPU | — | `qtt_ladder.json`; Re 100 re-run `results/rerun_Re100_20260914.json` |
| Dense vs tensor-train wall-clock and error, Re 10 → 10⁴ | CPU, workstation | — | `airbus_bench_workstation.json` |
| Spectral storage floor on the canonical flow | CPU | — | `floor_airbus_spectral.json` |

### Appendix B — Why the linear operator is exact for this benchmark

For u = U_c − cos(x − U_c t) sin y·e^{−2νt}, v = sin(x − U_c t) cos y·e^{−2νt}, the convective term (u·∇)u equals −∇(¼e^{−4νt}(cos 2(x − U_c t) + cos 2y)) + U_c ∂_x u, a gradient plus the mean-flow advection; the gradient is absorbed by the pressure, and what remains is u_t + U_c u_x = ν∇²u for each component. The finite-volume discretisation of that operator with central fluxes is the matrix A used by every solver in this package, so the Schrödingerised, dense and tensor-train routes are graded on one operator against one analytic solution.

### Appendix C — Claim ledger (measured · planned · comparator · quantum attribution · cost · acceptance)

| claim | status | classical comparator | quantum attribution | total cost charged | acceptance threshold |
|---|---|---|---|---|---|
| Schrödingerised solver reproduces the analytic vortex to finite-volume discretisation error | measured (exact simulation, 15–19 qubits + p register) | classical solve of the same operator; analytic solution | algorithm output, simulated | minutes CPU | error within 2× of FV at every rung |
| Memory 2·log₂N + n_p qubits vs 4N²·8 bytes | exact | grid storage | register size | — | — |
| Time-to-solution crossover at N ≈ 65,536 (Re 10⁴) | modelled, constants declared | explicit FV, identical stencils; pseudo-spectral named as SOTA for this flow | Hamiltonian simulation | readout charged at O(1/ε) | crossover < 2²⁰ per axis at ≤ 50 logical qubits under compiled constants (Phase 2) |
| 42-qubit evolution primitive on hardware | measured | exact statevector of the medium | hardware execution | 1 banked job | correlation ≥ 0.9 after one envelope (met) |
| Tensor-train annex: Re ladder, memory ratio vs dense | measured | dense FV; spectral floor (128 bytes) | quantum-inspired | workstation hours | L2 ≤ 10⁻⁴ per rung (met) |
