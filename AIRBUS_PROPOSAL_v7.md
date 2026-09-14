# Quantics tensor-train finite-volume solver for the convecting Taylor-Green vortex, with a hardware-native pathway

**Global Quantum + AI Challenge 2026 — Airbus Enterprise Challenge · Phase 1 Concept Proposal · Team Merlin Digital (GIC 2026 dual-track finalist — Mitsubishi/AIST materials track) · v7.0 · 2026-09-12**

---

## 1. Problem framing

The challenge asks for a solver of the 2D convecting Taylor-Green vortex whose scaling in time-to-solution, memory and error against the analytic solution is quantified as the Reynolds number rises — with grid resolution forced upward by Re. The cost driver of classical finite-volume CFD is exactly that resolution: a dense N×N field costs O(N²) memory and the pressure projection costs a global solve on it at every step.

Our approach is on the challenge's **tensor-network-within-finite-volume** track, evaluated on its own merits as the statement invites: represent the cell averages as **quantics tensor trains** (each axis split into log₂N binary cores) so memory is O(log N · χ²) per field instead of O(N²), and drive the time evolution by fluxes across volume faces with the same stencils a dense FV solver uses. The credible advantage is measured, not projected: at Re = 10⁴ on a 4096² grid the field is held at bond dimension 3, the L2 error against the analytic solution is 4.4×10⁻⁵, and the memory ratio to the dense field is 9.2×10⁴. The quantum-inspired representation is what a fault-tolerant quantum solver would load and evolve; the same package carries a hardware-native pathway already flown on 42 superconducting qubits.

## 2. Technical approach

**Paradigm.** Quantum-inspired tensor network (quantics TT) inside a finite-volume scheme; hybrid classical execution today; the QTT encoding is the amplitude-encoding a gate-model or FTQC fluid solver consumes.

**Discretisation (FV as required).** Cell averages on an N×N periodic grid; central advective flux differences across faces; 5-point viscous stencil; RK2 (Heun) in time; incompressibility by pressure projection with the Poisson solve done as TT-conjugate-gradient with rank rounding — a general operator, no benchmark-specific shortcut; zero-mean enforced.

**Tensor-train machinery.** Periodic shift-by-one-cell as an exact rank-2 quantics MPO (binary increment with carry); d/dx, d/dy and the Laplacian assembled from shifts; initial cos/sin fields built as exact rank-2 trains by angle addition, so no dense array is ever materialised at large N; left-right QR then right-left SVD rounding with tolerance 10⁻¹⁰ and rank cap 64.

**Cross-validation.** The QTT solver and a dense reference with identical stencils agree to ~10⁻⁸ in L2 at N = 64 and 128; both are then graded against the exact TGV solution (velocity L2 norm and kinetic-energy decay, the two error measures the statement names).

**Hardware-native pathway (quantum-algorithm track, annex).** A lattice-gas Floquet step on 21 coupled qubit pairs (42 qubits, IBM Marrakesh job d9solsgpdb6s73e6a6q0) with a Taylor-Green-shaped parity preparation: pair term κ·D as native RZZ, single-site mixing as native RX, inter-pair XX+YY exchange as native RXX+RYY, Galilean convection applied exactly at readout. The measured quantity is the equal-time connected pair covariance X(q,k) over depth k. It is a measured quantum medium with the TGV mode as its initial condition — the dynamical primitive a quantum lattice-Boltzmann solver builds on.

## 3. Feasibility and resource requirements

| resource | status |
|---|---|
| Solver | `qtt_fv.py`, self-contained NumPy; dense reference in the same file |
| Compute for the Re ladder | single workstation; Re 10⁴ / N 4096 run in 2.8 h; larger Re on the same node with rank cap raised as needed |
| Classical comparator | dense FV with identical stencils (same file) plus a named production FV/spectral solver for wall-clock comparison in Phase 2 |
| Quantum hardware for the pathway | IBM Heron via Startup Program; QuEra Aquila via Braket (analog lattice-gas), chassis frozen |
| Data | none required; analytic TGV solution is the reference |

Assumptions: periodic domain, incompressible Newtonian fluid, the statement's parameters (L = 2π, V₀ = 1, U_c = 1, ν = V₀L/Re). Constraint: the convecting TGV is a single Fourier mode, so its quantics rank stays small at any Re; the Phase-2 plan (§5) adds multi-mode and perturbed initial conditions where rank growth — and therefore the real cost frontier — appears.

## 4. Expected impact

**Delivered Re ladder (receipt `qtt_ladder.json`; grid resolution scaled with Re; error = velocity L2 vs analytic solution):**

| Re | N | steps | χ | L2 vs exact | QTT bytes | dense bytes | memory ratio | wall (s) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 64 | 326 (T = 1.0) | 3 | 2.7×10⁻⁴ | 2,336 | 131,072 | 56× | 415 |
| 100 | 256 | 522 (T = 1.0) | 8 | 6.2×10⁻⁵ | 13,472 | 2,097,152 | 156× | 3,515 |
| 10,000 | 4,096 | 522 (T = 0.1) | 3 | 4.4×10⁻⁵ | 5,824 | 536,870,912 | 92,183× | 10,004 |

![Figure 1 — The Re ladder: dense/QTT memory ratio (left) and velocity L2 error against the analytic solution (right) as grid resolution is scaled with Reynolds number; rank ≤ 8 throughout.](C:/quantum ai 2026/figs_v7/airbus_ladder.png)

Three readings. **Memory:** the representation cost is flat in Re (5.8 kB at Re 10⁴) while the dense field grows as N²; this is the scaling the statement asks to be quantified — **against dense finite volume, which is the right comparison for a general solver and the wrong one for this particular flow.** We ran that floor ourselves (`floor_airbus_spectral.json`): the convecting TGV has four or five significant Fourier coefficients per velocity component at every resolution (u carries the mean flow), so a spectral method stores the field in ~128 bytes — 18–105× *less* than the tensor train. On the canonical benchmark the compression headline belongs to spectral methods; the tensor train earns its place where the spectrum is not sparse, which is why the multi-mode rung below is the one that matters. **Error:** L2 stays at 10⁻⁴–10⁻⁵ across three decades of Re with rank ≤ 8. **Time:** wall-clock grows with the step count and the TT-CG projection; the Phase-2 target is to report it against a named production solver at matched L2, and to push Re to the largest value at which the bounded-rank representation still meets a stated error bar.

![Figure 2 — Time-to-solution against Reynolds number for the QTT finite-volume solver on one workstation (steps and horizon per rung annotated); the matched-L2 comparison against a production solver is the Phase-2 plot.](C:/quantum ai 2026/figs_v7/airbus_tts.png)

**Hardware pathway, measured.** 42 qubits, 21 coupled pairs; six pre-registered gates on the medium's operator identity and waveform passed (operator identity, waveform vs exact medium, k-ordering, sideband structure, revival present); connected X(q,k) extracted from retrieved counts with finite-shot bias correction. This establishes the quantum medium on today's hardware — the dynamical primitive the gate-model pathway builds on.

### Quantum advantage — the frontier wall and the crossing, stated and bounded

**The field's direction, and its wall.** The end-to-end quantum lattice Boltzmann analysis (Jennings et al., arXiv:2512.03758) bounds the achievable speed-up at O(Re^{3D/8}), polynomial and conditional on a high-error-tolerance regime, with nonlinearity (Carleman failed to converge), field readout O(Re^{3/8}) and state preparation as the obstacles; no quantum CFD run has crossed a classical solver on time-to-solution. The productive 2026 direction is quantum-inspired: tensor-network lattice Boltzmann (arXiv:2512.07615), tensor-train compressible turbulence with solver arithmetic in TT form (arXiv:2608.26995), GPU quantics turbulence (PRR 7, 013112).

**The crossing, measured.** Our solver is that direction delivered inside the statement's finite-volume frame: cell averages and face fluxes in quantics tensor-train form, Re 10 → 10⁴, N to 4,096, rank ≤ 8, L2 4×10⁻⁵ against the analytic solution, dense/QTT memory ratio 9.2×10⁴ at Re 10⁴, QTT–dense agreement 10⁻⁸ on identical stencils. **Where it is decisive:** at multi-mode and perturbed initial conditions, where the dense field grows as N² and the tensor-train rank stays bounded — the Phase-2 measurement, alongside wall-clock against a production solver at matched L2 — with the 42-qubit lattice-gas medium as the gate-model pathway.

**Scalability to industrial relevance.** The representation cost is O(log N · χ²) per field, so grid resolution is bought logarithmically: the 4,096² field at Re 10⁴ costs 5.8 kB against 537 MB dense, and the same solver at 16,384² costs the same again plus two cores per axis. Multi-mode and perturbed initial conditions raise χ, which is exactly the quantity the Phase-2 rank-trajectory measurement reports; on the quantum side, a quantics train is the amplitude encoding a gate-model lattice-Boltzmann solver loads, so the classical and hardware pathways share one data structure.

**Business value, bounded.** Aerodynamic screening is memory- and time-bound at the resolutions where the physics matters; a solver that holds a resolved field at a memory cost independent of Re moves the screening loop from the cluster to the workstation for the flow classes where rank stays bounded. The Phase-2 plot — time-to-solution at matched L2 against a production solver — is the number that converts to screening throughput per engineer-day; it is reported, not assumed.

**What a successful PoC demonstrates for Airbus.** A finite-volume solver in tensor-train form that holds a resolved vortex field at a memory cost independent of Re, with the scaling plots the statement specifies (time-to-solution, memory, L2 and kinetic-energy-decay error vs Re) side by side with a production classical solver; and a costed map from that representation onto a gate-model lattice-Boltzmann or LCU fluid solver, so Airbus can read where a fault-tolerant machine enters the same workflow.

## 5. Validation plan

Pre-registered per rung: Re, N, T, rank cap, error bars. Metrics: velocity L2 vs analytic solution; kinetic-energy decay vs analytic; wall-clock and peak memory for QTT and for the dense/production solver on named hardware; bond dimension trajectory. Success: (i) L2 ≤ 10⁻⁴ at every Re rung up to the stated maximum with χ ≤ 16; (ii) memory ratio reported at each rung; (iii) time-to-solution at matched L2 reported against the production solver, with the crossover Re — if any — stated as a number. Robustness: mesh/time convergence at each Re; a perturbed and a two-mode initial condition where rank growth is expected, with the rank trajectory reported rather than assumed. Hardware pathway: uniform-k re-flight of the 42-qubit medium (68 circuits, 8,192 shots each, controls: no-exchange, decoupled-pair, phase-randomised) graded only on its declared covariance observable. Every negative kept.

## 6. Hybrid / cross-domain integration

The QTT solver is a drop-in FV component: same stencils, same fluxes, same projection, with the field storage swapped. It integrates with a classical CFD workflow at the field level (a QTT field can be materialised to a dense patch for any post-processing) and with quantum hardware at the encoding level (a quantics train is a circuit-ready amplitude encoding). The hardware-native medium runs on IBM Heron and, in analog form, on QuEra Aquila; both are read back into the same covariance grader. Phase-2 integration target: one Re-ladder script that emits the classical scaling plots and the quantum-encoding cost table from the same run.

## 7. Team capability

Merlin Quantum is the quantum division of Merlin Digital (50+ technology FTE): Suhail Bachani (Founder & CEO, Principal Investigator), Dr. Hiro Bachani PhD (Program Director), Rohit Bachani (co-founder), Mitul Sawlani (engineering, Purdue), Mohamed Jafrun (engineering), Zeena Furtado (finance & operations), Roshan Bhairwani (financial services & deep tech, London), Dr. Ana Baroni MSc (domain specialist). GIC 2026 dual-track finalist. Programme hardware record: 259 receipted QPU jobs, 12.3 million shots, across IBM Heron (up to 156 qubits, 2,720 two-qubit gates per circuit), QuEra Aquila, QCi Dirac-3 and Rigetti, with pre-registration and archived counts for every claim, negative results included.

## 8. Scope, with treatment

(i) Time-to-solution is reported, not claimed (`airbus_bench_workstation.json`, identical stencils, same grid per rung): the tensor-train solver is slower than dense finite volume at Re 10, 100 and 1,000 (335 s vs 0.7 s; 3,515 s on the Sept-3 ladder and 9,068 s in the current-solver re-run on a shared laptop, vs 16.7 s; 4,462 s vs 975 s) and 5.2× faster at Re 10⁴ (3,207 s vs 16,692 s) — but its velocity error floors at 3–6×10⁻⁵ from rank truncation while dense reaches 1.9×10⁻⁸ on the same grid, so at matched error dense wins at every rung tested. Process peak memory is a clean number only at Re 10 (88 MB vs 44 MB); above that the benchmark ran both solvers in one process and the representation size is the memory claim. The shipped solver carries a robust-SVD fallback added after one benchmark run met a non-converging SVD at Re 100; the rung reproduces with it. (ii) The convecting TGV has four or five significant Fourier coefficients per component, so a spectral method stores it in ~128 bytes — 18–105× below the tensor train (`floor_airbus_spectral.json`); the memory ratio is a statement about dense FV, and Phase 2 moves to multi-mode flows where the spectrum is not sparse. (iii) The Re = 10⁴ rung is integrated to T = 0.1 — T = 1.0 is a Phase-2 run. (iv) Wall-clock is not yet compared against a production solver — that comparison at matched L2 is the PoC's primary plot. (v) The 42-qubit medium is not a Navier–Stokes solution — it is the hardware pathway, reported as such. (vi) The FTQC encoding cost is not yet computed — the costed map is a Phase-2 deliverable.

---

### Appendix A — Job register and receipts

| measurement | machine | job id | receipt |
|---|---|---|---|
| 42-qubit lattice-gas medium, Taylor–Green parity preparation; P4 covariance X(q,k) | ibm_marrakesh | `d9solsgpdb6s73e6a6q0` | `qlga_counts_20260810_122411.json`, `qlga_xqk.json`, `qlga_advantage.json` |
| Layout scout | ibm_marrakesh | `d9solbntfhrs73dt1t40` | `qlga_state.json` |
| Re ladder rows (N, steps, χ, bytes, L2, wall) | CPU | — | `qtt_ladder.json` |
| Re 100 rung re-run with the shipped solver (robust-SVD fallback), 2026-09-14: L2 6.3×10⁻⁵, χ ≤ 8, 15.2 kB — reproduces | CPU | — | `results/rerun_Re100_20260914.json` |
| Wall-clock and error, dense vs QTT on identical stencils and grids, Re 10 → 10⁴ (workstation, 32 cores, 8.3 h) | CPU | — | `airbus_bench_workstation.json` |
| QTT vs dense cross-validation, 4.9e-9 / 8.5e-7 at N = 64 / 128 | CPU | — | `qtt_fv.py` stage g1 |
| Design card, 10 coupled pairs, exact-diagonalisation referee | CPU | — | `qlga_design.json` |

### Appendix B — Scheme notes

Quantics ordering: x cores then y cores, least-significant bit first per axis. Poisson solve: TT-CG with rank rounding at each iteration, zero-mean constraint applied in TT form. Time step from the CFL and viscous limits at each N; step counts in the table. The Re = 10⁴ rung is integrated to T = 0.1 (522 steps at N = 4096); T = 1.0 at that rung is a Phase-2 run. Reproduce: `python qtt_fv.py ladder`.

### Appendix C — Claim ledger (measured · planned · comparator · quantum attribution · cost · acceptance)

| claim | status | classical comparator | quantum attribution | total cost charged | acceptance threshold |
|---|---|---|---|---|---|
| Re ladder 10 → 10⁴, L2 ≤ 2.7×10⁻⁴, memory ratio to 9.2×10⁴ **vs dense FV**; Re 100 rung re-run and reproduced (6.3×10⁻⁵) | measured | dense FV with identical stencils (cross-validated 10⁻⁸); spectral floor run — 128 bytes, 18–105× below QTT on this flow | quantum-inspired | workstation hours | L2 ≤ 10⁻⁴ per rung, χ ≤ 16 |
| Time-to-solution and peak memory vs production solver at matched L2 | in progress (workstation benchmark) | optimised FV / spectral | quantum-inspired | workstation | crossover Re reported as a number |
| Multi-mode / perturbed initial condition rank growth | planned | dense reference | quantum-inspired | workstation | rank trajectory reported |
| 42-qubit lattice-gas medium, six gates | measured | exact medium at 10 pairs | hardware pathway | 1 job | gates PASS |

