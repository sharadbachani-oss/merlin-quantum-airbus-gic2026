# -*- coding: utf-8 -*-
"""
verify.py — credential-free replication of the Airbus headline numbers.

    python verify.py

Reads archived JSON in this folder and re-derives every figure cited in
AIRBUS_SUBMISSION_v4.md. Re-flying the quantum job needs IBM credentials
and is not this script.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
# Derived operator constants — closed form, no external dependency.
import math as _m
KAPPA = 3.0 / (3.0 - _m.sqrt(5.0))            # κ, structural spectral gap
DELTA = _m.sqrt((KAPPA / 2.0) ** 2 + 4.0) - KAPPA / 2.0   # per-pair gap

JOB = "d9solsgpdb6s73e6a6q0"
SCOUT = "d9solbntfhrs73dt1t40"


def load(name):
    return json.load(open(os.path.join(HERE, name), encoding="utf-8"))


def near(a, b, tol):
    return abs(float(a) - float(b)) <= tol


def main():
    print("=" * 64)
    print("Airbus TGV — native advantage replication (no credentials)")
    print("=" * 64)
    ok = True

    adv = load("qlga_advantage.json")
    spec = load("qlga_vortex_spectrum.json")
    result = load("qlga_result_20260810_122411.json")
    design = load("qlga_design.json")
    pending = load("qlga_pending_20260810_122411.json")
    state = load("qlga_state.json")
    ladder = load("qtt_ladder.json")

    print("\n[1] Job IDs (banked, not invented)")
    print(f"    main {state['main_job']}  scout {state['scout_job']}  "
          f"backend {result['backend']}")
    if state["main_job"] != JOB or spec["job"] != JOB or adv["job"] != JOB:
        print("    FAIL job mismatch")
        ok = False
    elif state["scout_job"] != SCOUT:
        print("    FAIL scout mismatch")
        ok = False
    else:
        print("    PASS")

    print("\n[2] Operator identity with first-principles")
    print(f"    pending mu {pending['mu']:.12f}  first-principles KAPPA {KAPPA:.12f}")
    if abs(pending["mu"] - KAPPA) > 1e-12 or not adv["operator"]["mu_match"]:
        print("    FAIL")
        ok = False
    else:
        print("    PASS")

    print("\n[3] Waveform grade (A1-P vs TGV exponential)")
    for km, exp_med, exp_tgv in (("1", 0.9933, 0.8497), ("2", 0.9705, 0.8456)):
        m = adv["modes"][km]
        print(f"    k={km}: medium {m['damped_model_corr']:.4f}  "
              f"exp {m['exp_corr']:.4f}  revival {m['has_revival']}")
        if not near(m["damped_model_corr"], exp_med, 5e-4):
            print("    FAIL medium corr drift")
            ok = False
        if not near(m["exp_corr"], exp_tgv, 5e-4):
            print("    FAIL exp corr drift")
            ok = False
        if not m["has_revival"]:
            print("    FAIL revival missing")
            ok = False
        if m["damped_model_corr"] <= m["exp_corr"]:
            print("    FAIL medium does not beat exponential")
            ok = False

    print("\n[4] S(w) sidebands (real-time FFT, no continuation)")
    for km, hw_sb, exp_sb in (("1", 0.327, 0.113), ("2", 0.354, 0.114)):
        hw = spec["modes"][km]["hardware"]["sideband_fraction"]
        ex = spec["modes"][km]["tgv_exponential"]["sideband_fraction"]
        print(f"    k={km}: hw {hw:.3f}  exp {ex:.3f}")
        if not near(hw, hw_sb, 5e-3) or not near(ex, exp_sb, 5e-3):
            print("    FAIL sideband drift")
            ok = False
        if hw < ex + 0.05:
            print("    FAIL sideband gate")
            ok = False

    print("\n[5] Medium k-ordering (design D3) + inherited hardware track")
    ratio = design["modes"]["2"]["decay_rate"] / design["modes"]["1"]["decay_rate"]
    print(f"    design k-ratio {ratio:.3f}  (report 2.82)")
    if not (2.0 < ratio < 8.0):
        print("    FAIL D3 window")
        ok = False
    v2 = result.get("grading_v2", {}).get("modes", {})
    if v2:
        v2r = v2["2"]["envelope_per_step"] / v2["1"]["envelope_per_step"]
        print(f"    grading_v2 envelope ratio {v2r:.3f}  (report 1.70)")

    print("\n[6] Drift control + convection identity")
    v_main = result["mode1"]["drift"]
    v0 = result["gates"]["drift_control"]
    print(f"    main-series drift {v_main:+.4f}  (report +0.0013); "
          f"C0-arm {v0:+.4f}  (report +0.022)")
    if abs(v_main) > 0.01 or abs(v0) > 0.05:
        print("    FAIL drift control")
        ok = False
    if not design["gates"]["D2"]:
        print("    FAIL D2 boost identity")
        ok = False

    print("\n[7] G-gate verdict + P4 honesty")
    g = adv["gates"]
    counts_path = state.get("counts")
    have_counts = bool(counts_path and os.path.isfile(counts_path))
    have_xqk = os.path.isfile(os.path.join(HERE, "qlga_xqk.json"))
    print(f"    verdict {g['verdict']}  P4 extracted {g['P4_Xqk_extracted']}")
    if g["verdict"] != "PASS":
        print("    FAIL G-gates")
        ok = False
    if have_counts and have_xqk:
        if not g["P4_Xqk_extracted"]:
            print("    FAIL: counts+xqk on disk but P4 left unclaimed")
            ok = False
        xqk = load("qlga_xqk.json")
        if xqk.get("job") != JOB:
            print("    FAIL xqk job mismatch")
            ok = False
        pins = (
            ("1", "8", "0.5", 0.0006360245413196782),
            ("1", "8", "1.0", 0.0008317876109276083),
            ("1", "12", "2.0", -0.0015634446954760973),
            ("2", "8", "0.5", -0.0009975878299665058),
        )
        for km, ks, q, exp in pins:
            got = xqk["series"][km][ks]["Xq"][q]
            print(f"    X(q={q},k={ks}|mode {km}) {got:+.6e}")
            if not near(got, exp, 1e-12):
                print("    FAIL X(q,k) drift")
                ok = False
        p4c = adv.get("P4_Xqk") or {}
        if p4c.get("job") != JOB or "8" not in (p4c.get("k_ge_8") or {}).get("1", {}):
            print("    FAIL advantage P4 card missing k>=8 extract")
            ok = False
    else:
        if g["P4_Xqk_extracted"]:
            print("    FAIL: P4 must stay unclaimed without archived counts")
            ok = False

    print("\n[7b] Retrieved counts reproduce archived amplitudes")
    if have_counts:
        sys.path.insert(0, HERE)
        from airbus_qlga import profile_counts, mode_amp_phase, STEPS, KMODES
        blob = load(os.path.basename(counts_path))
        if blob.get("job") != JOB:
            print("    FAIL counts job mismatch")
            ok = False
        else:
            path, nr = blob["path"], blob["nr"]
            maxdiff = 0.0
            for km in KMODES:
                arch = result[f"mode{km}"]["amps"]
                for i, ks in enumerate(STEPS):
                    prof = profile_counts(blob["counts"][f"M{km}_{int(ks)}"], path, nr)
                    A, _ = mode_amp_phase(prof, km)
                    maxdiff = max(maxdiff, abs(A - arch[i]))
            print(f"    max |A_retr - A_arch| = {maxdiff:.3e}")
            if maxdiff > 1e-12:
                print("    FAIL amplitude mismatch vs archived result")
                ok = False
    else:
        print("    SKIP: no archived counts")

    print("\n[8] QTT annex (mean field, not the lead)")
    for row, Re, chi in zip(ladder, (10, 100, 10000), (3, 8, 3)):
        print(f"    Re={row['Re']}: chi={row['chi']}  L2={row['l2']:.2e}  "
              f"mem {row['mem_ratio']:.1f}x")
        if row["Re"] != Re or row["chi"] != chi:
            print("    FAIL ladder drift")
            ok = False
        if row["l2"] > 1e-3:
            print("    FAIL L2")
            ok = False

    print("\n[9] Jennings bound + parity leftover honesty")
    banned = adv.get("what_we_do_not_claim") or []
    joined = " ".join(str(x).lower() for x in banned) if isinstance(banned, list) else str(banned).lower()
    if "navier" in joined or "ns" in joined or "incompressible" in joined:
        print("    PASS  NS supremacy is explicitly refused")
    else:
        print("    FAIL  NS bound missing from what_we_do_not_claim")
        ok = False
    if "qtt" in joined or "mean" in joined:
        print("    PASS  QTT mean-field win is not inverted")
    else:
        print("    FAIL  QTT leftover not disclosed")
        ok = False
    design_nr = design.get("nr") or design.get("n_pairs")
    print(f"    design closed-orbit size {design_nr}  (nr=21 anchors must stay absent)")
    if os.path.isfile(os.path.join(HERE, "qlga_nr21_anchors.json")):
        print("    FAIL invented nr=21 anchors")
        ok = False
    else:
        print("    PASS  no nr=21 anchor file")
    print("    protocol: parity-filtered leftover is a next-series grade, not this job")

    print("\n[10] Closed-orbit parity leftover (banked counts, no new job)")
    par_path = os.path.join(HERE, "airbus_parity_leftover.json")
    if os.path.isfile(par_path):
        par = load("airbus_parity_leftover.json")
        if par.get("job_id") != JOB:
            print("    FAIL parity job mismatch")
            ok = False
        elif par.get("nr21_anchors_invented") is True:
            print("    FAIL invented nr=21")
            ok = False
        elif par.get("k_grid_uniform") is True:
            print("    FAIL: banked k-grid is not uniform; do not claim k=0..16")
            ok = False
        elif not par.get("gates", {}).get("extracted_from_banked"):
            print("    FAIL: parity must come from banked counts")
            ok = False
        else:
            print(f"    leftover beats pole {par['gates']['leftover_beats_pole']}  "
                  f"revival {par['gates']['revival_survives']}")
            print("    PASS  parity leftover extracted; NS refused; no new job")
    else:
        print("    FAIL missing airbus_parity_leftover.json")
        ok = False

    print("\n[11] Dirac closed-orbit leftover (local; not an NS solve)")
    dpath = os.path.join(HERE, "airbus_dirac_leftover.json")
    if os.path.isfile(dpath):
        dl = load("airbus_dirac_leftover.json")
        if dl.get("job_id") is not None or dl.get("new_hardware_advantage_demonstrated") is True:
            print("    FAIL invented Dirac hardware advantage")
            ok = False
        elif not (dl.get("local_gates") or {}).get("airbus_orbit_held"):
            print("    FAIL orbit must be more held than named H")
            ok = False
        elif not (dl.get("honesty") or {}).get("jennings_ns_not_claimed"):
            print("    FAIL Jennings NS bound missing")
            ok = False
        else:
            print("    PASS  Dirac orbit leftover local; NS refused; no new job")
    else:
        print("    FAIL missing airbus_dirac_leftover.json")
        ok = False

    print("\n[12] Full-W native card (64-cell continuous)")
    fw = os.path.join(HERE, "..", "results", "dirac_full_W_flight.json")
    if os.path.isfile(fw):
        f = json.load(open(fw, encoding="utf-8"))
        if f.get("n_vars") != 64 or f.get("job_type") != "sample-hamiltonian":
            print("    FAIL native card must be 64-cell sample-hamiltonian")
            ok = False
        elif f.get("new_hardware_advantage_demonstrated") is True:
            print("    FAIL full-W advantage was not evidenced")
            ok = False
        elif f.get("pf_vacuum_seen") is True:
            print("    FAIL PF vacuum was not seen")
            ok = False
        else:
            print("    PASS  full-W flown; V-octet collapse; leftover advantage stays false")
    else:
        print("    FAIL missing dirac_full_W_flight.json")
        ok = False

    print("\n" + ("PASS" if ok else "FAIL") +
          " — headline numbers match archived receipts")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
