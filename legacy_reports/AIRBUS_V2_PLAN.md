# Airbus Submission v2 — Work Plan (2026-08-10)

**Base: the existing GPQLBM package (SUBMISSION.pdf, April 2026).
Verdict on v1: correct kernel, honest claims matrix, good structure — but
simulator-only (no hardware receipt anywhere) and mitigation-benefit
framed as the headline instead of advantage. v2 fixes both, using the
program's now-validated machinery. Same portal limits (5 files).**

## The two upgrades that change the submission's class

### Upgrade 1 — GPQLBM ON REAL HARDWARE (IBM, free open plan)
v1's gauge_pass = 0.86 and the entire noise model are inherited constants;
the "ideal stabilizer detection" caveat is the submission's weakest line.
Fix by flying the smallest honest instance:
- N=2..3 grid of 6-qubit GPQLBM cells (24–54 qubits) on ibm_fez/kingston,
  1–3 Trotter steps, Z-basis readout.
- MEASURE: (a) real per-cell gauge pass rates on TGV workload (replaces
  the inherited 0.86), (b) real detection efficiency (inject known
  single-qubit errors on flagged circuits — the eraser/einselection
  machinery), (c) one hardware KE point vs analytic TGV with error bars.
- Flight discipline as always: frozen pending, scout->fly->grade,
  model==artifact statevector gate at 24q, receipts published (not NDA).
- Result: "first hardware-executed gauge-diagnosed QLBM timestep, with
  measured — not assumed — detection performance." No other entrant will
  have receipts.

### Upgrade 2 — TN-in-FV SOLVER (the challenge's second track) with
### MEASURED quantum-inspired advantage at high Re
The challenge explicitly invites Tensor Networks within a Finite Volume
framework and counts quantum-inspired advantage TODAY. This is where a
real, measured advantage claim lives (the week's lesson: claim advantage
where the crossover is measurable, not projected):
- Solver: cell-AVERAGE fields (per track spec) as an MPS over the x-index
  (y folded, or MPO-2D), flux-based FV update (Rusanov/central flux for
  advection + viscous flux), pressure projection via FFT-diagonal
  Poisson in TT form; bond dimension chi adaptive.
- The advantage measurement: memory and time-to-solution vs a dense FV
  reference at matched L2 error, swept Re = 10 -> 100 -> 1000 -> 10^4
  with Reynolds-required N. The smooth TGV field compresses at low chi;
  chi(Re) growth IS the scaling story, measured not projected — and the
  box's validated TEBD/MPS engines make this a fast build.
- Deliverables per the challenge: time-to-solution, memory, L2/KE error
  scaling curves vs Re; dense-FV classical baseline on identical grids.

## Keep from v1 (already good)
- D2Q9 kernel + analytic benchmark harness (validated to 0.3–0.4% KE).
- Claims/withhold discipline and structure.
- Aquila geometry as the analog annex (fly only if budget approved —
  costs real money, unlike IBM).

## Order of work
1. TN-FV solver core + dense FV baseline (box GPU/CPU; engine exists).
2. IBM GPQLBM mini-card: emit + statevector-gate the 24q instance (here).
3. Fly the mini-card (open plan, 2 jobs: scout + main).
4. Re-scaling sweep on TN-FV (box); classical baseline table.
5. Rewrite SUBMISSION v2: advantage exhibit first (TN-FV measured
   crossover + hardware receipts), mitigation story second, claims matrix
   updated with measured constants.

## Deadline
Same program as VW (Phase I closes Sept 15). Airbus and VW submissions
share the September window; both are writing-complete by Sept 8.
