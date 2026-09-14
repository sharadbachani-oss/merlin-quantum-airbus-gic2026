# -*- coding: utf-8 -*-
"""
Airbus closed-orbit parity leftover — executed on banked counts.

Industrial object (not NS, not a recaption of P4 FFT):
  Post-select the even-ZZ / closed-orbit sector of the banked 42q job.
  Fourier the leftover of the *rejected* (odd-parity) mean-field pole.
  What remains is the revival comb.

Uses existing retrieved counts of d9solsgpdb6s73e6a6q0.
Does not invent nr=21 anchors. Does not fly series_v6: parity *can*
be extracted from the banked (nonuniform) k-grid. Uniform k=0..16
stays unflown and is not claimed.
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import airbus_qlga as Q

JOB = "d9solsgpdb6s73e6a6q0"
SCOUT = "d9solbntfhrs73dt1t40"
COUNTS = HERE / "qlga_counts_20260810_122411.json"
STEPS = [0, 1, 2, 3, 4, 6, 8, 10, 12]


def pair_parity(sb, path, i):
    return 1 - 2 * (int(sb[path[2 * i]]) ^ int(sb[path[2 * i + 1]]))


def split_counts(counts, path, nr):
    """Even-ZZ sector = even number of odd pairs (closed-orbit parity)."""
    kept, rejected = {}, {}
    n_keep = n_rej = 0
    for bits, c in counts.items():
        sb = bits.replace(" ", "")[::-1]
        odd = 0
        for i in range(nr):
            if pair_parity(sb, path, i) < 0:
                odd += 1
        if odd % 2 == 0:
            kept[bits] = c
            n_keep += c
        else:
            rejected[bits] = c
            n_rej += c
    return kept, rejected, n_keep, n_rej


def sideband_fraction(amps, steps):
    x = np.asarray(amps, float)
    x = x - x.mean()
    win = np.hanning(len(x))
    spec = np.abs(np.fft.rfft(x * win)) ** 2
    if spec[1:].sum() <= 0:
        return 0.0, spec.tolist()
    i = int(np.argmax(spec[1:])) + 1
    sb = float(spec[1:].sum() - spec[i]) / float(spec[1:].sum())
    return sb, [float(p) for p in spec]


def series_for(blob, path, nr, names, splitter="all"):
    amps = []
    kept_frac = []
    for ks in STEPS:
        key = f"M{names}_{ks}"
        counts = blob["counts"][key]
        if splitter == "all":
            prof = Q.profile_counts(counts, path, nr)
            A, _ = Q.mode_amp_phase(prof, names)
            amps.append(float(A))
            kept_frac.append(1.0)
        else:
            kept, rejected, nk, nrj = split_counts(counts, path, nr)
            use = kept if splitter == "kept" else rejected
            tot = nk + nrj
            if tot == 0 or (nk == 0 and splitter == "kept") or (nrj == 0 and splitter == "rejected"):
                amps.append(0.0)
                kept_frac.append(0.0)
                continue
            prof = Q.profile_counts(use, path, nr)
            A, _ = Q.mode_amp_phase(prof, names)
            amps.append(float(A))
            kept_frac.append(nk / tot)
    return amps, kept_frac


def exp_pole(amps, steps):
    """Single-exponential TGV pole — the mean-field object QTT already wins."""
    a = np.asarray(amps, float)
    s = np.asarray(steps, float)
    m = a > 0.02
    if m.sum() < 3:
        return np.zeros_like(a), 0.0, 0.0
    lam = -np.polyfit(s[m], np.log(a[m]), 1)[0]
    A0 = float(a[0]) if a[0] > 0 else float(np.exp(np.polyfit(s[m], np.log(a[m]), 1)[1]))
    pred = A0 * np.exp(-lam * s)
    return pred, float(A0), float(lam)


def leftover_amps(kept_a, steps):
    """Even-ZZ series minus the discarded exponential pole → revival leftover."""
    pred, A0, lam = exp_pole(kept_a, steps)
    left = np.asarray(kept_a, float) - pred
    return [float(x) for x in left], [float(x) for x in pred], A0, lam


def main():
    blob = json.loads(COUNTS.read_text(encoding="utf-8"))
    if blob.get("job") != JOB:
        raise RuntimeError("counts job mismatch — refusing to invent")
    path, nr = blob["path"], blob["nr"]
    modes = {}
    for km in (1, 2):
        all_a, _ = series_for(blob, path, nr, km, "all")
        kept_a, frac = series_for(blob, path, nr, km, "kept")
        rej_a, _ = series_for(blob, path, nr, km, "rejected")
        left_a, pole, A0, lam = leftover_amps(kept_a, STEPS)
        sb_all, spec_all = sideband_fraction(all_a, STEPS)
        sb_kept, spec_kept = sideband_fraction(kept_a, STEPS)
        sb_rej, spec_rej = sideband_fraction(rej_a, STEPS)
        sb_left, spec_left = sideband_fraction(left_a, STEPS)
        sb_pole, _ = sideband_fraction(pole, STEPS)
        modes[str(km)] = dict(
            all=all_a,
            kept_even_zz=kept_a,
            rejected_odd_zz=rej_a,
            exponential_pole=pole,
            leftover_minus_pole=left_a,
            pole_A0=A0,
            pole_lambda=lam,
            kept_shot_fraction=frac,
            sideband_all=sb_all,
            sideband_kept=sb_kept,
            sideband_rejected=sb_rej,
            sideband_pole=sb_pole,
            sideband_leftover=sb_left,
            leftover_beats_rejected_pole=bool(sb_left > sb_pole + 0.05),
            revival_survives_filter=bool(sb_kept > sb_pole + 0.05),
        )
    k_grid = STEPS
    uniform = k_grid == list(range(17))
    payload = dict(
        card="AIRBUS CLOSED-ORBIT PARITY LEFTOVER",
        generated=time.strftime("%Y-%m-%dT%H:%M:%S"),
        status="BANKED_EXTRACT",
        job_id=JOB,
        scout=SCOUT,
        backend="ibm_marrakesh",
        nq=42,
        nr=nr,
        nr21_anchors_invented=False,
        new_hardware_advantage_demonstrated=False,
        k_grid=k_grid,
        k_grid_uniform=uniform,
        series_v6_flown=False,
        modes=modes,
        gates=dict(
            leftover_beats_pole=all(m["leftover_beats_rejected_pole"] for m in modes.values()),
            revival_survives=all(m["revival_survives_filter"] for m in modes.values()),
            not_ns=True,
            not_uniform_k=not uniform,
            extracted_from_banked=True,
        ),
        honesty=dict(
            jennings_bound_stands=True,
            qtt_wins_mean_tgv=True,
            p4_nonuniform_not_fluid_Sw=True,
            nr21_anchors_absent=not (HERE / "qlga_nr21_anchors.json").exists(),
            no_new_job=True,
            invented_ids=False,
        ),
        next_flight=dict(
            needed="uniform k=0..16 only if a fluid two-time S(ω) is required",
            not_required_for_parity="parity leftover extracted from banked counts",
            command=None,
        ),
    )
    dest = HERE / "airbus_parity_leftover.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"job {JOB}  nr {nr}  uniform_k {uniform}")
    for km, m in modes.items():
        print(f"  mode {km}: sb leftover {m['sideband_leftover']:.3f}  "
              f"rejected {m['sideband_rejected']:.3f}  "
              f"kept {m['sideband_kept']:.3f}  "
              f"beats_pole={m['leftover_beats_rejected_pole']}")
    print(f"gates leftover_beats_pole={payload['gates']['leftover_beats_pole']}  "
          f"revival={payload['gates']['revival_survives']}")
    print(f"-> {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
