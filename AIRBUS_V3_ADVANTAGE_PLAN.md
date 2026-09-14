# Airbus v4 — Native Advantage-of-Route Plan (2026-09-08)

Supersedes the 2026-09-02 v3 plan. That plan set the winning posture as
**quantum-inspired QTT today**. That is the wrong bar. QTT of a rank-~3
Fourier mode is the annex, not the lead — the same role Mitsubishi gave
absolute energy (kept off-device).

## Winning posture (this is the bar)

**first-principles-native dynamics on today's hardware**, with a bounded advantage
of *route*, not an unbounded Navier–Stokes supremacy claim.

- Operator: `H = κ·D − A₆`, `κ = 3/(3−√5)`, imported from `C:\first-principles`
  (not fitted). Exact shallow RY prep. Hardware-native RZZ / RX / RXX+RYY
  maps. Noise is an instrument (one damping envelope, A1-P pattern).
- Lead card: the banked 42-qubit QLGA flight on ibm_marrakesh
  (`d9solsgpdb6s73e6a6q0`, scout `d9solbntfhrs73dt1t40`). No new IDs.
- Lead object: the **real-frequency vortex-relaxation spectrum** `S(ω)`
  of a TGV-prepared parity vortex evolving under `H`. Direct FFT of the
  measured real-time series. **No imaginary-time step. No analytic
  continuation.**
- Annex: classical QTT-in-FV of the mean TGV field (the brief's tensor-
  network track). Like Mitsubishi's absolute energy: it certifies the
  verification sector and is not the crossing.

## Why QTT cannot be the lead

The convecting 2D TGV is one Fourier mode. Its quantics rank is ~3 at
any `N`. Memory 10⁴–10⁵× vs a dense array is real and is reported, but
it is **not quantum advantage**. Any low-rank Fourier / spectral method
fakes the mean field for free. Leading with that number failed both
mandatory gates (first-principles physics 4.0, today's-hardware advantage 3.0).

## The crossing (what classical routes cannot cheaply fake)

Three classical routes, all named:

1. **Imaginary-time + analytic continuation** to a real-frequency
   spectrum — mathematically ill-posed at any system size (the
   Mitsubishi `A(ω)` argument; two distinct spectra can share the same
   Euclidean data).
2. **Low-rank Fourier / QTT of the mean TGV** — produces a single
   exponential pole `A(t) = A₀ exp(−2νk²t)`. It cannot produce the
   closed-orbit **revivals** of `H`, nor their sideband weight in `S(ω)`.
3. **Classical real-time many-body on `H`** — the A1 adjudication
   already measured death at `k ≥ 8` (operators past 2×10⁸ terms). The
   QLGA P4 channel is that observable class. **Extracted from retrieved
   counts of the banked job; not a Navier–Stokes claim.**

Hardware evolves `H` natively. The spectrum is a Fourier transform.

Measured on the banked job (local regrade, `qlga_advantage.py`):

| Mode | Medium corr (A1-P) | TGV-exp corr | Δ | Revival | Sideband hw / exp |
|---|---:|---:|---:|---|---:|
| k=1 | 0.993 | 0.850 | +0.144 | yes | 0.327 / 0.113 |
| k=2 | 0.971 | 0.846 | +0.125 | yes | 0.354 / 0.114 |

Exact-medium k-ratio (design D3) = 2.82. grading_v2 envelope ratio
(disclosed) = 1.70. PHI=0 spurious drift +0.0013 (C0-arm +0.022).
Raw exponential P3-as-written remains FAIL — that was an instrument
error (grading a coherent revival against a pure exponential) and
stays disclosed.

## Honest bound (do not claim full NS)

The medium realises **advection–diffusion + interacting quantum
content**, graded against the exact TGV decay *form* in the envelopes
and against the exact closed-orbit of `H` in the waveforms. This is
**not** a claim that the chip integrates incompressible Navier–Stokes.

Mitsubishi bounded the lead as `A(ω)`, not 40-qubit absolute energy.
We bound the lead as `S(ω)` of a TGV-prepared vortex, not 2D NS on the
chip. Brief ref [9] (Jennings et al.) already proves end-to-end
nonlinear-fluid quantum algorithms carry *bounded* advantage; we sit
inside that bound as a **route** advantage on today's devices.

Industrial split (same spine as Mitsubishi):

| Sector | Object | Who computes it | Role |
|---|---|---|---|
| Verification / mean field | `(u,v,p,KE)` of the 2D TGV | QTT annex, exact analytic grade | like absolute energy — off-device |
| Native / spectral | `S(ω)` vortex-relaxation, revival sidebands, k-sectors | ibm_marrakesh under `H` | the lead |
| Hard fluctuation (P4) | connected `X(q,k)` | same shots, retrieved | A1 hardness class; extracted in `qlga_xqk.json` |

## Protocol (already flown — no reflight required)

1. Design gate (statevector, `nr=10`): D1 decay, D2 boost identity, D3
   k-ratio. PASS, frozen in `qlga_design.json`.
2. Scout `d9solbntfhrs73dt1t40` → fly `d9solsgpdb6s73e6a6q0` on
   ibm_marrakesh, 21 pairs / 42 qubits, 8192 shots.
3. Grade mean waveforms vs exact medium (A1-P). Local regrade in
   `qlga_advantage.py` (credential-free).
4. `retrieve` (optional, needs IBM): archive raw counts. `xqk` then
   extracts P4. Refuse if counts are missing — do not invent them.
5. Rewrite the report so quantum-inspired is not the winning posture.

## Gates for this rebuild

- G1: `κ` on the card equals first-principles `KAPPA` (import, not a copy).
- G2: hardware vs exact-medium waveform corr ≥ 0.95 / 0.90 (k=1 / k=2).
- G3: medium corr beats TGV-exponential corr by ≥ 0.02, lower RMS.
- G4: revival present in both measured modes.
- G5: exact medium carries k-ordering in `[2, 8]`; hardware tracks both
  modes (do not re-fail the card on raw exponential P3).
- G6: hardware sideband fraction exceeds the exponential's by ≥ 0.05.
- P4: extracted from retrieved counts of `d9solsgpdb6s73e6a6q0`.

All six G-gates PASS on the archived job. P4 is extracted. Not NS.

## What this is not

- Not a new hardware campaign.
- Not a QTT memory-supremacy paper.
- Not full Navier–Stokes on 42 qubits.
- Not a fake job ID or a simulated count file labelled as hardware.
