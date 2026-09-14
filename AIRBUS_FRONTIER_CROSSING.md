# Airbus — closed-orbit parity leftover (frontier crossing)

**11 September 2026.** Novel industrial object: the **enstrophy-cascade leftover after a closed-orbit / ring-parity filter** — revival sidebands that survive the filter, not a recaption of vortex \(S(\omega)\) and not incompressible NS.

Status: **REAL on the banked medium** + **BANKED EXTRACT** of the parity leftover (`airbus_parity_leftover.json`) + **LOCAL DIRAC** closed-orbit \(L(\tau)\). Uniform \(k=0..16\) still unflown. Jennings NS bound stays.
`new_hardware_advantage_demonstrated: false` (v6 contract). Banked job is prior. No new job.

---

## SOTA wall (what currently wins)

| Winner today | What it produces | What it cannot cheaply produce |
|---|---|---|
| QTT / tensor-train in finite volume (Gourianov *et al.*; this package `qtt_fv.py`) | Mean TGV to \(10^{-9}\)–\(10^{-7}\). Ladder: Re=10 \(\chi=3\) L2 \(2.70\times10^{-4}\); Re=100 \(\chi=8\) L2 \(6.18\times10^{-5}\); Re=10\(^4\) \(\chi=3\) L2 \(4.37\times10^{-5}\) | Closed-orbit **revivals** and \(S(\omega)\) sidebands of the interacting medium |
| Low-rank Fourier / single-pole TGV exponential | A decaying mean mode | Sideband fraction 0.327 / 0.354 (hardware) vs 0.113 / 0.114 (exponential) |
| DNS / LES / LBM; Carleman–LB FTQC proposals | \((u,v,p)\) on a mesh | A NISQ claim of incompressible NS. **Jennings et al. 2025** (arXiv:2512.03758, brief [9]) **bounds** end-to-end fluid QA: modest, high error-tolerance, \(\mathcal{O}(\mathrm{Re}^{3/4(1+D/2)}\,q_M)\) with \(q_M=\mathcal{O}(\mathrm{Re}^{3/8})\) for drag; numerical D=2 closer to \(\mathcal{O}(\mathrm{Re}^{1.936}q_M)\) |
| Imaginary-time + analytic continuation of a vortex correlator | A smooth pole | The realtime revival comb |

QTT is the correct attack on the **mean** TGV and it **wins** there. That is why it is the annex. Using “QTT cannot see revivals” as same-task supremacy is a category error.

---

## Novel crossing — parity-filtered leftover, not another FFT of \(X(q,k)\)

Banked lead (keep, do not re-caption): ibm_marrakesh `d9solsgpdb6s73e6a6q0` (42q / 21 pairs, 8192 shots) + scout `d9solbntfhrs73dt1t40`. Hardware vs exact-medium corr **0.993 / 0.971**; vs TGV exponential **0.850 / 0.846**. Operator \(H=\mu^\star D-A_6\) plus inter-pair XX+YY. Exact shallow prep. Floquet step, not automatically \(\exp(-iHt)\).

The **new** object, inside the Jennings cap:

**Closed-orbit parity leftover.** The medium’s conserved ring / even-ZZ parity is a noise filter. Post-select (or syndrome-weight) the closed-orbit sector; Fourier the leftover of the *rejected* mean-field pole. What remains is the revival comb — the enstrophy-cascade spectral content of the lattice gas. Low-rank Fourier of the mean TGV is exactly the discarded pole. Continuation of \(G(\tau)\) cannot put the sidebands back.

**Wrap-before-cone** (secondary, protocol): onset-of-turbulence as the depth where the vortex light-cone wraps the box. Do not claim that wrap is NS instability.

P4 connected \(X(q,k)\) from the same shots (mode-1 \(k=8\): \(X(0.5)=6.36\times10^{-4}\), \(X(1.0)=8.32\times10^{-4}\)) inherits A1 \(k\ge 8\) death at \(2\times10^8\) terms. v6: retrieved \(k\)-grid is **nonuniform**; ordinary FFT of P4 is **rejected** as a fluid two-time structure factor. Uniform \(k=0..16\) is a separate unflown series (`series_v6`, 68 circuits).

---

## Why today’s hardware can bound a *route*

- 21-pair heavy-hex map, exact prep, one banked flight already reproduces the closed-orbit theorem after one damping envelope (`nr=10` design; **`nr=21` box anchors are not in this folder**).
- Noise damping is the envelope we fit, not something we ZNE into a mean-field win.
- Jennings already says the NS mesh is not the target. We stay on a 42-qubit spectral / fluctuation object.

---

## Honesty bound

From `qlga_advantage.json` `what_we_do_not_claim`:

- Full incompressible NS on the chip.
- Mean-field L2 / KE supremacy over QTT (QTT wins).
- Unbounded fluid QA.
- That `nr=21` anchors exist here.
- That P4 on the nonuniform grid is a fluid \(S(\omega)\).
- Simulated counts as hardware. No new job IDs.

---

## Dirac remap (11 Sep) — closed-orbit leftover, not an NS solve

Dirac card is allowed only as a **new leftover / parity object**. Even-parity occupation is the closed-orbit filter; \(L(\tau)\) of the hex orbit graph is the leftover of the rejected mean-field pole. Jennings 2025 still bounds end-to-end fluid QA. Do not caption this as incompressible NS.

**Executed 11 Sep (local 2⁶ Gibbs):** closed-orbit \(L(1)=0.0164\) (most held of the six graphs); even-parity filter holds at \(\tau=4\). `airbus_dirac_leftover.json`.

**Shared Dirac card (flown, crossing not evidenced):** `6aa32614…`–`…6f30`. Jennings NS bound stays. Track `job_id: null`. Advantage **false**.

## Next flight (protocol only)

Shared Dirac leftover card only if the readout is leftover/parity, not a claimed NS solve. Uniform \(k=0..16\) IBM series stays unflown.

**Executed 10 Sep:** parity leftover extracted from banked counts of `d9solsgpdb6s73e6a6q0` (`python airbus_parity_leftover.py`). Even-ZZ kept series minus the single-exponential TGV pole. Mode-1 leftover sideband 0.450 vs pole; mode-2 0.666 vs pole; leftover-beats-pole **true**; revival-survives-filter **true**. k-grid remains nonuniform — **not** claimed as uniform \(k=0..16\). nr=21 anchors **not invented**. No new flight (parity *was* extractable).

Local gate: `python verify.py` `[9]` Jennings honesty, `[10]` parity leftover, `[11]` Dirac orbit leftover.
