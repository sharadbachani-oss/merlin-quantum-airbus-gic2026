# Airbus TGV — Merlin Digital (GIC 2026)

Native quantum lattice-gas evolution of a TGV-prepared vortex on today's
IBM Heron hardware, plus a classical QTT-in-FV annex of the mean field.

**Read first:** `AIRBUS_SUBMISSION_v4.md` (live text after P4 retrieve). `AIRBUS_SUBMISSION_v3.md` is the pre-retrieve snapshot.
**Advantage claim:** `AIRBUS_ADVANTAGE_EXHIBIT.md`.
**Plan:** `AIRBUS_V3_ADVANTAGE_PLAN.md` (v4 posture; QTT is not the lead).

## Lead result

ibm_marrakesh `d9solsgpdb6s73e6a6q0` (42 qubits). Operator
`H = κ·D − A₆` imported from `C:\our first-principles engine`. Deliverable is the
vortex-relaxation spectrum `S(ω)` — real-time FFT, no analytic
continuation. QTT memory compression of the rank-~3 TGV is the annex.

## Verify (no credentials)

```
pip install numpy
python verify.py
```

Expected: G1–G6 PASS, P4 extracted from retrieved counts, job IDs match the banked flight.

## Layout

| File | Role |
|---|---|
| `airbus_qlga.py` | native protocol (design/scout/fly/grade/retrieve/xqk/advantage) |
| `qlga_advantage.py` | local spectrum + route grade |
| `qtt_fv.py` | TN-in-FV annex |
| `qlga_result_20260810_122411.json` | archived hardware amplitudes |
| `qlga_counts_20260810_122411.json` | retrieved raw counts (banked job) |
| `qlga_xqk.json` | P4 \(X(q,k)\) extract |
| `qlga_vortex_spectrum.json` | `S(ω)` artifact |
| `RECEIPTS.md` | claim → file → job ID |

## Bound

Advection–diffusion + interacting medium. Not incompressible
Navier–Stokes on the chip. Brief ref [9] already bounds end-to-end
nonlinear-fluid quantum advantage; this is a route advantage inside
that bound.
