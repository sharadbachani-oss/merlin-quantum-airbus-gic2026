"""Airbus: the operator a fluid actually has, and the comparator on it.

THE DIAGNOSIS THAT LED HERE. The previous formulation put the coupled pairs on a
one-dimensional nearest-neighbour chain. A matrix-product state follows that at
bond dimension 4, flat in system size, so no advantage exists on it. That was a
fault in the formulation, not in the physics: a fluid perturbation field is not
a chain. It is a field on a two-dimensional surface, and shock-boundary-layer
interaction couples spanwise as well as chordwise.

THE FIX. Place each coupled pair on a site of a two-dimensional grid and put
exchange bonds on both lattice directions. The local generator is unchanged,
kappa*D - sum X per pair with structural kappa. Only the connectivity changes,
from a line to a surface.

WHY THIS IS THE DECIDING VARIABLE. A matrix-product state has to thread a
one-dimensional path through a two-dimensional lattice, so any cut crosses the
transverse width and the bond dimension carries that width. This is the standard
reason matrix-product methods are strong in one dimension and weak in two. The
test below is controlled: identical qubit count, chain against grid, so the only
difference is connectivity.

SCOPE. Measured cost for these instances, this observable, this accuracy target
and the orderings searched. A credible classical baseline, not a complexity
proof.
"""
import json, math, os
import numpy as np

KAPPA = 3.0 / (3.0 - math.sqrt(5.0))
GB = 1.0
FIDELITY = 0.999
HERE = os.path.dirname(os.path.abspath(__file__))


def grid_bonds(Lx, Ly):
    """Pair p at grid (x,y) occupies qubits 2p, 2p+1. Exchange bonds join
    neighbouring pairs in BOTH lattice directions."""
    site = lambda x, y: y * Lx + x
    bonds = []
    for y in range(Ly):
        for x in range(Lx):
            if x + 1 < Lx:
                bonds.append((site(x, y), site(x + 1, y)))
            if y + 1 < Ly:
                bonds.append((site(x, y), site(x, y + 1)))
    return bonds


def chain_bonds(npairs):
    return [(i, i + 1) for i in range(npairs - 1)]


def build(npairs, bonds, g=GB):
    """diag part and the off-diagonal action, matrix-free."""
    n = 2 * npairs
    D = 1 << n
    idx = np.arange(D, dtype=np.int64)
    diag = np.zeros(D)
    for p in range(npairs):
        za = 1 - 2 * ((idx >> (2 * p)) & 1)
        zbb = 1 - 2 * ((idx >> (2 * p + 1)) & 1)
        diag += 0.5 * (1.0 - za * zbb)
    diag *= KAPPA
    # exchange acts between the second qubit of one pair and the first of the next
    xb = [(2 * a + 1, 2 * b) for a, b in bonds]
    return n, diag, xb, g


def apply_H(psi, n, diag, xb, g):
    out = diag * psi
    v = psi.reshape([2] * n)
    for i in range(n):
        out -= np.flip(v, axis=i).reshape(-1)
    idx = np.arange(psi.size, dtype=np.int64)
    for a, b in xb:
        fa, fb = 1 << a, 1 << b
        unlike = ((idx >> a) & 1) != ((idx >> b) & 1)
        src = idx[unlike]
        out[src ^ fa ^ fb] += 2.0 * g * psi[src]
    return out


def ground(n, diag, xb, g, iters=300):
    rng = np.random.default_rng(21)
    v = rng.normal(size=1 << n); v /= np.linalg.norm(v)
    shift = KAPPA * n + n + 4.0 * g * len(xb) + 4.0
    for _ in range(iters):
        v = shift * v - apply_H(v, n, diag, xb, g)
        v /= np.linalg.norm(v)
    return v


def lanczos(psi, n, diag, xb, g, dt, m=16):
    V = np.empty((m, psi.size), dtype=np.complex128)
    a = np.zeros(m); b = np.zeros(m)
    beta = np.linalg.norm(psi); V[0] = psi / beta; j = m
    for i in range(m):
        w = apply_H(V[i], n, diag, xb, g)
        a[i] = np.vdot(V[i], w).real
        Vi = V[: i + 1]
        w -= Vi.T @ (Vi.conj() @ w)
        nb = np.linalg.norm(w)
        if nb < 1e-13:
            j = i + 1; break
        b[i] = nb
        if i + 1 < m: V[i + 1] = w / nb
    T = np.diag(a[:j]) + np.diag(b[: j - 1], 1) + np.diag(b[: j - 1], -1)
    ev, U = np.linalg.eigh(T)
    return beta * ((U @ (np.exp(-1j * ev * dt) * U[0].conj())) @ V[:j])


def min_bond(psi, n, order):
    v = psi.reshape([2] * n)
    v = np.transpose(v, [n - 1 - q for q in order])
    need = 1
    for cut in range(1, n):
        s = np.linalg.svd(v.reshape(1 << cut, -1), compute_uv=False)
        p = s ** 2; p = p / p.sum()
        need = max(need, min(int(np.searchsorted(np.cumsum(p), FIDELITY) + 1), p.size))
    return need


def orderings(npairs, Lx=None, Ly=None):
    n = 2 * npairs
    o = {"sequence": list(range(n))}
    pa = []
    for p in range(npairs):
        pa += [2 * p, 2 * p + 1]
    o["pair_adjacent"] = pa
    if Lx and Ly:                       # snake the grid, the standard MPS choice
        snake = []
        for y in range(Ly):
            xs = range(Lx) if y % 2 == 0 else reversed(range(Lx))
            for x in xs:
                p = y * Lx + x
                snake += [2 * p, 2 * p + 1]
        o["snake"] = snake
    return o


def measure(label, npairs, bonds, Lx=None, Ly=None, T=6.0, npts=7):
    n, diag_nat, _, g = build(npairs, [])
    n, diag, xb, g = build(npairs, bonds)
    psi = ground(n, diag_nat, [], g).astype(np.complex128)
    dt = T / (npts - 1)
    best = {k: 1 for k in orderings(npairs, Lx, Ly)}
    for _ in range(npts - 1):
        psi = lanczos(psi, n, diag, xb, g, dt)
        for k, o in orderings(npairs, Lx, Ly).items():
            best[k] = max(best[k], min_bond(psi, n, o))
    cap = 1 << (n // 2)
    chi = min(best.values())
    qcx = 8 * len(bonds) + npairs
    print("  %-14s pairs=%2d sites=%2d bonds=%2d  chi=%s  best=%3d of cap %4d  qCX=%d"
          % (label, npairs, n, len(bonds), best, chi, cap, qcx), flush=True)
    return dict(label=label, npairs=npairs, n_sites=n, n_bonds=len(bonds),
                chi_by_ordering=best, chi_best=chi, chi_cap=cap,
                quantum_logical_cx=qcx)


if __name__ == "__main__":
    print("controlled test: identical qubit count, chain against grid")
    rows = []
    for (Lx, Ly) in [(2, 2), (3, 2), (4, 2)]:
        npairs = Lx * Ly
        rows.append(measure("grid %dx%d" % (Lx, Ly), npairs, grid_bonds(Lx, Ly), Lx, Ly))
        rows.append(measure("chain %d" % npairs, npairs, chain_bonds(npairs)))
    print()
    for (Lx, Ly) in [(2, 2), (3, 2), (4, 2)]:
        npairs = Lx * Ly
        gr = [r for r in rows if r["label"] == "grid %dx%d" % (Lx, Ly)][0]
        ch = [r for r in rows if r["label"] == "chain %d" % npairs][0]
        print("  %2d pairs: grid chi=%3d  chain chi=%3d   grid/chain = %.2fx"
              % (npairs, gr["chi_best"], ch["chi_best"],
                 gr["chi_best"] / ch["chi_best"]))
    grids = [r for r in rows if r["label"].startswith("grid")]
    chains = [r for r in rows if r["label"].startswith("chain")]
    gg = [grids[i]["chi_best"] / grids[i - 1]["chi_best"] for i in range(1, len(grids))]
    cc = [chains[i]["chi_best"] / chains[i - 1]["chi_best"] for i in range(1, len(chains))]
    print("\n  chi growth per added pair-column:  grid %s   chain %s"
          % (["%.2fx" % x for x in gg], ["%.2fx" % x for x in cc]))
    json.dump(dict(card="fluid operator: two-dimensional connectivity vs chain",
                   fidelity_target=FIDELITY, rows=rows,
                   grid_chi_growth=gg, chain_chi_growth=cc,
                   scope="measured cost for these instances, observable, accuracy "
                         "target and orderings searched; credible classical "
                         "baseline, not a complexity proof"),
              open(os.path.join(HERE, "results", "grid_operator_comparator.json"), "w"),
              indent=1)
    print("  wrote results/grid_operator_comparator.json")
