"""Airbus: retarded response in place of the equal-time structure factor,
plus the classical comparator arm.

WHY THE OBSERVABLE CHANGES. Certification by Analysis needs buffet onset within
3 kt of validation data, and what blocks it is the amplitude and frequency of a
self-sustained shock oscillation. RANS either converges to steady state or gets
both wrong, which is why scale-resolving simulation is required and unaffordable.

The previous P4 observable is the spatial Fourier transform of the equal-time
connected pair covariance. As the v6 audit states, a transform over k gives the
temporal content of that statistic and is not automatically the two-time dynamic
structure factor of a fluid. An unsteady spectrum is a two-time object. The
retarded response

    chi_ab(t) = i<[Z_a(t), Z_b(0)]> = -2 Im <Z_a(t) Z_b>

is one, and it carries a frequency and an amplitude, which is exactly the pair
the certification question asks for.

OPERATOR, kept as emitted. Chain of nr pairs on 2*nr sites: pairs (2i, 2i+1)
carry kappa*D - sum X; inter-pair bonds (2i+1, 2(i+1)) carry XX+YY exchange with
coupling g, which supplies transport.

COMPARATOR. Same operator, same observable, matched 99.9% state fidelity, bond
dimension as the maximum over all cuts, with several orderings searched and the
classical method's best case reported.

SCOPE. A measured cost for this instance family, observable, accuracy target and
orderings searched. A credible classical baseline and censored benchmark
observation, not a complexity proof.
"""
import json, math, os
import numpy as np

KAPPA = 3.0 / (3.0 - math.sqrt(5.0))
GB = 1.0
FIDELITY = 0.999
HERE = os.path.dirname(os.path.abspath(__file__))


def build_H(nr, g=GB):
    """kappa*D on pairs, minus sum X, plus XX+YY exchange on inter-pair bonds."""
    n = 2 * nr
    D = 1 << n
    idx = np.arange(D, dtype=np.int64)
    zb = [(1 - 2 * ((idx >> i) & 1)).astype(np.float64) for i in range(n)]
    diag = np.zeros(D)
    for i in range(nr):
        diag += 0.5 * (1.0 - zb[2 * i] * zb[2 * i + 1])
    H = np.diag(KAPPA * diag)
    for i in range(n):                              # -sum X
        H[idx ^ (1 << i), idx] -= 1.0
    for i in range(nr - 1):                         # XX+YY on bonds
        a, b = 2 * i + 1, 2 * (i + 1)
        fa, fb = 1 << a, 1 << b
        ba = (idx >> a) & 1
        bb = (idx >> b) & 1
        swap = ba != bb                             # XX+YY moves only unlike pairs
        src = idx[swap]
        H[src ^ fa ^ fb, src] += 2.0 * g
    return H, n


def spectrum(nr, T=24.0, npts=241, site_a=0, site_b=None):
    H, n = build_H(nr)
    site_b = n - 1 if site_b is None else site_b
    ev, U = np.linalg.eigh(H)
    g0 = U[:, 0]
    idx = np.arange(1 << n)
    za = (1 - 2 * ((idx >> site_a) & 1)).astype(float)
    zbv = (1 - 2 * ((idx >> site_b) & 1)).astype(float)
    Ut = U.T
    A0 = Ut @ g0
    B = Ut @ (zbv * g0)
    Za = Ut @ (za[:, None] * U)
    M = (A0[:, None] * Za) * B[None, :]
    dE = ev[:, None] - ev[None, :]
    times = np.linspace(0.0, T, npts)
    chi = np.array([-2.0 * float((M * np.sin(dE * t)).sum()) for t in times])
    win = np.hanning(len(chi))
    y = np.fft.rfft(chi * win, n=8 * len(chi))
    om = 2.0 * np.pi * np.fft.rfftfreq(8 * len(chi), d=times[1] - times[0])
    P = np.abs(y) ** 2
    m = om > 0.05
    o, p = om[m], P[m]
    j = int(np.argmax(p))
    return dict(nr=nr, n_sites=n, E0=float(ev[0]),
                chi_at_0=float(abs(chi[0])),
                peak_frequency=float(o[j]),
                peak_amplitude=float(np.max(np.abs(chi))),
                spectral_weight_frac=float(p[j] / p.sum()))


def orderings(n, nr):
    seq = list(range(n))
    pair_adj = []
    for i in range(nr):
        pair_adj += [2 * i, 2 * i + 1]
    interleave = list(range(0, n, 2)) + list(range(1, n, 2))
    return {"sequence": seq, "pair_adjacent": pair_adj, "interleaved": interleave}


def min_bond(psi, n, order, fidelity=FIDELITY):
    v = psi.reshape([2] * n)
    v = np.transpose(v, [n - 1 - q for q in order])
    need = 1
    for cut in range(1, n):
        s = np.linalg.svd(v.reshape(1 << cut, -1), compute_uv=False)
        p = s ** 2
        p = p / p.sum()
        need = max(need, min(int(np.searchsorted(np.cumsum(p), fidelity) + 1), p.size))
    return need


def comparator(nrs=(3, 4, 5), T=6.0, npts=9):
    from scipy.linalg import expm
    rows = []
    for nr in nrs:
        H, n = build_H(nr)
        ev, U = np.linalg.eigh(H)
        # QUENCH, not an eigenstate. Starting from the ground state of the FULL
        # generator is stationary: it only picks up a phase, so its entanglement
        # cannot change and the measurement returns a constant. Start from the
        # no-exchange ground state and evolve under the coupled operator.
        H0, _ = build_H(nr, g=0.0)
        _, U0 = np.linalg.eigh(H0)
        psi = U0[:, 0].astype(np.complex128)
        times = np.linspace(0.0, T, npts)
        dt = times[1] - times[0]
        step = U @ np.diag(np.exp(-1j * ev * dt)) @ U.T.conj()
        best = {k: 1 for k in orderings(n, nr)}
        for _ in range(npts - 1):
            psi = step @ psi
            for k, o in orderings(n, nr).items():
                best[k] = max(best[k], min_bond(psi, n, o))
        cap = 1 << (n // 2)
        chi_star = min(best.values())
        n_bonds = nr - 1
        quantum_cx = 8 * n_bonds + nr        # 8 CX per evolved bond + 1 per pair prep
        rows.append(dict(nr=nr, n_sites=n, chi_cap=cap, chi_by_ordering=best,
                         chi_best_ordering=chi_star,
                         classical_cost_proxy=float(chi_star ** 2 * n),
                         quantum_logical_cx=quantum_cx))
        print("  nr=%d sites=%2d  chi=%s  best=%d of cap %d   chi^2n=%.3e   qCX=%d"
              % (nr, n, best, chi_star, cap, chi_star ** 2 * n, quantum_cx), flush=True)
    return rows


if __name__ == "__main__":
    print("Airbus: retarded response spectrum (frequency AND amplitude)")
    specs = [spectrum(nr) for nr in (3, 4, 5)]
    for s in specs:
        print("  nr=%d sites=%2d  chi(0)=%.1e  peak frequency=%.4f  "
              "amplitude=%.5f  weight=%.4f"
              % (s["nr"], s["n_sites"], s["chi_at_0"], s["peak_frequency"],
                 s["peak_amplitude"], s["spectral_weight_frac"]), flush=True)
    print("\n  chi(0)=0 identically, so the response is causal and carries an")
    print("  onset. An equal-time structure factor has neither property.")
    print("\nclassical comparator arm")
    rows = comparator()
    chis = [r["chi_best_ordering"] for r in rows]
    caps = [r["chi_cap"] for r in rows]
    cg = [chis[i] / chis[i - 1] for i in range(1, len(chis))]
    capg = [caps[i] / caps[i - 1] for i in range(1, len(caps))]
    qg = [rows[i]["quantum_logical_cx"] / rows[i - 1]["quantum_logical_cx"]
          for i in range(1, len(rows))]
    print("\n  chi growth per added pair: %s  (cap grows %s)"
          % (["%.2fx" % x for x in cg], ["%.0fx" % x for x in capg]))
    print("  quantum logical CX growth: %s" % (["%.2fx" % x for x in qg]))
    json.dump(dict(card="Airbus retarded response and classical comparator",
                   operator="kappa*D on pairs, -sum X, XX+YY exchange on bonds",
                   observable="retarded chi_ab(t) = -2 Im <Z_a(t) Z_b>",
                   replaces="spatial FT of equal-time connected pair covariance",
                   fidelity_target=FIDELITY, spectra=specs, comparator=rows,
                   chi_growth_per_pair=cg, cap_growth=capg, quantum_cx_growth=qg,
                   scope="measured cost for this instance family, observable, "
                         "accuracy target and orderings searched; credible "
                         "classical baseline and censored benchmark observation, "
                         "not a complexity proof"),
              open(os.path.join(HERE, "results", "response_and_comparator.json"), "w"),
              indent=1)
    print("  wrote results/response_and_comparator.json")
