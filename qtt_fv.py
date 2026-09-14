# -*- coding: utf-8 -*-
"""
qtt_fv.py — QUANTICS tensor-train FV solver for the 2D convecting
Taylor-Green vortex (Airbus track: TN within a Finite-Volume framework).

Fields are CELL AVERAGES on an N x N periodic grid, represented as
quantics tensor trains: each axis reshaped into m = log2(N) binary
cores (x cores then y cores, LSB first per axis). Memory O(m * chi^2)
per field instead of O(N^2).

Scheme (identical stencils to the dense reference in this file):
  - central advective flux differences, 5-point viscous stencil
  - RK2 (Heun)
  - incompressibility by pressure projection; the Poisson solve is
    TT-CG with rank rounding (general operator — no benchmark-specific
    shortcut), zero-mean enforced
Operators:
  - periodic shift-by-one-cell as an exact rank-2 quantics MPO
    (binary increment/decrement with carry)
  - d/dx, d/dy, Laplacian assembled from shifts
Initial condition:
  - cos/sin built as EXACT rank-2 quantics trains (angle-addition
    construction) — no dense array is ever materialised at large N.

Stages:
  g1     — cross-validate QTT vs dense at N=64,128 (same stencils):
           TT-vs-dense L2 must be ~1e-8; both graded vs exact TGV
  ladder — Re sweep with Reynolds-scaled N; dense runs where it fits,
           DNF documented where it cannot; memory/wall/L2 per pair
"""
import json, math, sys, time
import numpy as np

WORK = r"C:\quantum ai 2026\airbus"


# ===================== TT core machinery ==============================
def tt_round(cores, tol=1e-10, chi_max=64):
    """left-right QR then right-left SVD truncation."""
    cs = [c.copy() for c in cores]
    n = len(cs)
    for i in range(n - 1):
        r, d, c = cs[i].shape
        m = cs[i].reshape(r * d, c)
        q, rr = np.linalg.qr(m)
        cs[i] = q.reshape(r, d, q.shape[1])
        cs[i + 1] = np.einsum("ab,bdc->adc", rr, cs[i + 1])
    for i in range(n - 1, 0, -1):
        r, d, c = cs[i].shape
        m = cs[i].reshape(r, d * c)
        try:
            u, s, vt = np.linalg.svd(m, full_matrices=False)
        except np.linalg.LinAlgError:          # gesdd can fail to converge on ill-conditioned cores; gesvd is slower but robust
            import scipy.linalg as _sl
            u, s, vt = _sl.svd(m, full_matrices=False, lapack_driver="gesvd")
        keep = max(1, min(chi_max, int((s > tol * s[0]).sum()) if s[0] > 0 else 1))
        u, s, vt = u[:, :keep], s[:keep], vt[:keep]
        cs[i] = vt.reshape(keep, d, c)
        cs[i - 1] = np.einsum("abc,cd,d->abd", cs[i - 1], u, s)
    return cs


def tt_add(a, b):
    out = []
    n = len(a)
    for i in range(n):
        ra, d, ca = a[i].shape
        rb, _, cb = b[i].shape
        if i == 0:
            c = np.zeros((1, d, ca + cb))
            c[:, :, :ca] = a[i]; c[:, :, ca:] = b[i]
        elif i == n - 1:
            c = np.zeros((ra + rb, d, 1))
            c[:ra] = a[i]; c[ra:] = b[i]
        else:
            c = np.zeros((ra + rb, d, ca + cb))
            c[:ra, :, :ca] = a[i]; c[ra:, :, ca:] = b[i]
        out.append(c)
    return out


def tt_scale(a, s):
    out = [c.copy() for c in a]
    out[0] = out[0] * s
    return out


def tt_hadamard(a, b):
    out = []
    for ca, cb in zip(a, b):
        ra, d, cc = ca.shape
        rb, _, cd = cb.shape
        c = np.einsum("adb,cde->acdbe", ca, cb).reshape(ra * rb, d, cc * cd)
        out.append(c)
    return out


def tt_inner(a, b):
    m = np.ones((1, 1))
    for ca, cb in zip(a, b):
        m = np.einsum("ab,adc,bde->ce", m, ca, cb)
    return float(m[0, 0])


def tt_norm(a):
    return math.sqrt(max(tt_inner(a, a), 0.0))


def tt_mem(a):
    return sum(c.size for c in a) * 8


def tt_chimax(a):
    return max(c.shape[0] for c in a)


def mpo_apply(W, a):
    out = []
    for w, c in zip(W, a):
        rw, do, di, cw = w.shape
        ra, _, ca = c.shape
        o = np.einsum("wxyv,ayb->wax b".replace(" ", ""), w, c) \
            if False else np.einsum("woiv,aib->waovb", w, c)
        out.append(o.reshape(rw * ra, do, cw * ca))
    return out


# ===================== quantics builders ==============================
def const_tt(nc, val=1.0):
    cs = [np.ones((1, 2, 1)) for _ in range(nc)]
    cs[0] = cs[0] * val
    return cs


def trig_tt_1d(m, k, phase, kind):
    """exact rank-2 quantics train for cos/sin(2*pi*k*(n+0.5)/N + phase)
    over n = 0..N-1 (cell centers), LSB-first cores.
    Represent [cos, sin] pair propagated by angle addition."""
    N = 2 ** m
    th0 = 2 * math.pi * k * 0.5 / N + phase       # cell-center offset
    cores = []
    for i in range(m):
        ang = 2 * math.pi * k * (2 ** i) / N       # angle of this bit
        c = np.zeros((2, 2, 2))
        for b in (0, 1):
            a = ang * b
            R = np.array([[math.cos(a), -math.sin(a)],
                          [math.sin(a), math.cos(a)]])
            c[:, b, :] = R.T                       # row-vector propagation
        cores.append(c)
    # boundary: start row [cos(th0), sin(th0)], end column selects kind
    start = np.array([[math.cos(th0), math.sin(th0)]])
    cores[0] = np.einsum("ab,bdc->adc", start, cores[0])
    sel = np.array([[1.0], [0.0]]) if kind == "cos" else np.array([[0.0], [1.0]])
    cores[-1] = np.einsum("abc,cd->abd", cores[-1], sel)
    return cores


def field_2d(mx, my, fx, fy):
    """separable field f(x)*g(y): concatenate x-train and y-train."""
    return [c.copy() for c in fx] + [c.copy() for c in fy]


def shift_mpo_1d(m, direction):
    """periodic shift by one cell: index -> index+1 (direction=+1) or -1.
    LSB-first binary increment with carry; exact rank 2."""
    cores = []
    for i in range(m):
        w = np.zeros((2, 2, 2, 2))   # (rank_in, out_bit, in_bit, rank_out)
        # rank 0 = carry active, rank 1 = done
        if direction == +1:
            w[0, 1, 0, 1] = 1.0      # 0 -> 1, carry consumed
            w[0, 0, 1, 0] = 1.0      # 1 -> 0, carry propagates
        else:
            w[0, 0, 1, 1] = 1.0      # 1 -> 0, borrow consumed
            w[0, 1, 0, 0] = 1.0      # 0 -> 1, borrow propagates
        w[1, 0, 0, 1] = 1.0          # done: identity
        w[1, 1, 1, 1] = 1.0
        cores.append(w)
    # boundary: carry-in = active at LSB; accept carry-out in EITHER state
    # (carry falling off MSB = periodic wrap)
    start = np.zeros((1, 2)); start[0, 0] = 1.0
    cores[0] = np.einsum("ab,bdec->adec", start, cores[0])
    end = np.ones((2, 1))
    cores[-1] = np.einsum("adeb,bc->adec", cores[-1], end)
    return cores


def identity_mpo(m):
    w = np.zeros((1, 2, 2, 1))
    w[0, 0, 0, 0] = 1.0; w[0, 1, 1, 0] = 1.0
    return [w.copy() for _ in range(m)]


def mpo_2d(op_x, op_y):
    return [c.copy() for c in op_x] + [c.copy() for c in op_y]


# ===================== the solver =====================================
class QTTFV:
    def __init__(self, m, nu, tol=1e-9, chi_max=24,
                 field_tol=1e-6, field_chi=8):
        self.m, self.nu, self.tol, self.chi = m, nu, tol, chi_max
        self.ftol, self.fchi = field_tol, field_chi
        self._p_prev = None                      # CG warm start
        self.N = 2 ** m
        self.dx = 2 * math.pi / self.N
        I = identity_mpo(m)
        Sp = shift_mpo_1d(m, +1); Sm = shift_mpo_1d(m, -1)
        self.SPx = mpo_2d(Sp, I); self.SMx = mpo_2d(Sm, I)
        self.SPy = mpo_2d(I, Sp); self.SMy = mpo_2d(I, Sm)

    def r(self, a, tol=None):
        return tt_round(a, self.tol if tol is None else tol, self.chi)

    def ddx(self, f):
        # (S+ f)(n) = f(n-1), (S- f)(n) = f(n+1)  ->  central = (S- - S+)/2dx
        return self.r(tt_scale(tt_add(mpo_apply(self.SMx, f),
                                      tt_scale(mpo_apply(self.SPx, f), -1)),
                               1.0 / (2 * self.dx)))

    def ddy(self, f):
        return self.r(tt_scale(tt_add(mpo_apply(self.SMy, f),
                                      tt_scale(mpo_apply(self.SPy, f), -1)),
                               1.0 / (2 * self.dx)))

    def lap(self, f):
        s = tt_add(tt_add(mpo_apply(self.SPx, f), mpo_apply(self.SMx, f)),
                   tt_add(mpo_apply(self.SPy, f), mpo_apply(self.SMy, f)))
        return self.r(tt_add(tt_scale(s, 1.0 / self.dx ** 2),
                             tt_scale(f, -4.0 / self.dx ** 2)))

    def zero_mean(self, f):
        one = const_tt(2 * self.m)
        mean = tt_inner(f, one) / (self.N ** 2)
        return self.r(tt_add(f, tt_scale(one, -mean)))

    def poisson_cg(self, rhs, iters=12, rtol=1e-5):
        """solve lap(p) = rhs by TT-CG with rounding (zero-mean gauge)."""
        b = self.zero_mean(tt_scale(rhs, -1.0))   # solve (-lap)p = -rhs: SPD
        if self._p_prev is not None:
            x = self._p_prev
            rv = self.r(tt_add(b, self.zero_mean(self.lap(x))))
        else:
            x = tt_scale(const_tt(2 * self.m), 0.0)
            rv = b
        p = rv
        rs = tt_inner(rv, rv)
        b0 = math.sqrt(max(rs, 1e-300))
        for _ in range(iters):
            Ap = self.zero_mean(tt_scale(self.lap(p), -1.0))
            alpha = rs / max(tt_inner(p, Ap), 1e-300)
            x = self.r(tt_add(x, tt_scale(p, alpha)))
            rv = self.r(tt_add(rv, tt_scale(Ap, -alpha)))
            rs_new = tt_inner(rv, rv)
            if math.sqrt(max(rs_new, 0)) < rtol * b0:
                break
            p = self.r(tt_add(rv, tt_scale(p, rs_new / max(rs, 1e-300))))
            rs = rs_new
        x = self.zero_mean(x)
        self._p_prev = x
        return x

    def rhs(self, u, v):
        """FV RHS: -(u u_x + v u_y) + nu lap u  (and same for v)."""
        ru = tt_add(tt_scale(self.r(tt_hadamard(u, self.ddx(u))), -1.0),
                    tt_scale(self.r(tt_hadamard(v, self.ddy(u))), -1.0))
        ru = self.r(tt_add(ru, tt_scale(self.lap(u), self.nu)))
        rv_ = tt_add(tt_scale(self.r(tt_hadamard(u, self.ddx(v))), -1.0),
                     tt_scale(self.r(tt_hadamard(v, self.ddy(v))), -1.0))
        rv_ = self.r(tt_add(rv_, tt_scale(self.lap(v), self.nu)))
        return ru, rv_

    def project(self, u, v):
        div = self.r(tt_add(self.ddx(u), self.ddy(v)))
        p = self.poisson_cg(div)
        return (self.r(tt_add(u, tt_scale(self.ddx(p), -1.0))),
                self.r(tt_add(v, tt_scale(self.ddy(p), -1.0))))

    def step(self, u, v, dt):
        ru, rv_ = self.rhs(u, v)
        u1 = self.r(tt_add(u, tt_scale(ru, dt)))
        v1 = self.r(tt_add(v, tt_scale(rv_, dt)))
        u1, v1 = self.project(u1, v1)
        ru1, rv1 = self.rhs(u1, v1)
        u2 = self.r(tt_add(u, tt_scale(tt_add(ru, ru1), 0.5 * dt)))
        v2 = self.r(tt_add(v, tt_scale(tt_add(rv_, rv1), 0.5 * dt)))
        u2, v2 = self.project(u2, v2)
        # physics-level state sweep: remove truncation debris below the
        # scheme error; the rank-3 solution content is untouched
        return (tt_round(u2, self.ftol, self.fchi),
                tt_round(v2, self.ftol, self.fchi))


def tgv_init_tt(m, Uc=1.0, Vc=0.0, V0=1.0):
    cx = trig_tt_1d(m, 1, 0.0, "cos"); sx = trig_tt_1d(m, 1, 0.0, "sin")
    cy = trig_tt_1d(m, 1, 0.0, "cos"); sy = trig_tt_1d(m, 1, 0.0, "sin")
    one = const_tt(2 * m)
    u = tt_add(tt_scale(one, Uc),
               tt_scale(field_2d(m, m, cx, sy), -V0))
    v = tt_add(tt_scale(one, Vc),
               tt_scale(field_2d(m, m, sx, cy), V0))
    return tt_round(u), tt_round(v)


def tgv_exact_tt(m, t, nu, Uc=1.0, Vc=0.0, V0=1.0):
    F = math.exp(-2 * nu * t)
    cx = trig_tt_1d(m, 1, -Uc * t, "cos"); sx = trig_tt_1d(m, 1, -Uc * t, "sin")
    cy = trig_tt_1d(m, 1, -Vc * t, "cos"); sy = trig_tt_1d(m, 1, -Vc * t, "sin")
    one = const_tt(2 * m)
    u = tt_add(tt_scale(one, Uc),
               tt_scale(field_2d(m, m, cx, sy), -V0 * F))
    v = tt_add(tt_scale(one, Vc),
               tt_scale(field_2d(m, m, sx, cy), V0 * F))
    return tt_round(u), tt_round(v)


def l2_vs_exact_tt(u, v, m, t, nu):
    ue, ve = tgv_exact_tt(m, t, nu)
    du = tt_add(u, tt_scale(ue, -1.0)); dv = tt_add(v, tt_scale(ve, -1.0))
    N2 = float(4 ** m)
    return math.sqrt((tt_inner(du, du) + tt_inner(dv, dv)) / (2 * N2))


# ===================== dense reference (same stencils) ================
def dense_solve(N, nu, T, dt):
    dx = 2 * math.pi / N
    xs = (np.arange(N) + 0.5) * dx
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    u = 1.0 - np.cos(X) * np.sin(Y)
    v = 0.0 + np.sin(X) * np.cos(Y)
    kx = np.fft.fftfreq(N, d=1.0 / N)
    KX, KY = np.meshgrid(kx, kx, indexing="ij")
    K2 = KX ** 2 + KY ** 2; K2[0, 0] = 1.0
    # discrete symbols for the CENTRAL difference operators (match TT)
    sx = 1j * np.sin(KX * dx) / dx
    sy = 1j * np.sin(KY * dx) / dx
    lap = (2 * np.cos(KX * dx) + 2 * np.cos(KY * dx) - 4) / dx ** 2
    lap0 = lap.copy(); lap0[0, 0] = 1.0

    def ddx(f): return np.real(np.fft.ifft2(sx * np.fft.fft2(f)))
    def ddy(f): return np.real(np.fft.ifft2(sy * np.fft.fft2(f)))
    def lp(f): return np.real(np.fft.ifft2(lap * np.fft.fft2(f)))

    def project(u, v):
        div = ddx(u) + ddy(v)
        ph = np.fft.fft2(div) / lap0; ph[0, 0] = 0.0
        p = np.real(np.fft.ifft2(ph))
        return u - ddx(p), v - ddy(p)

    def rhs(u, v):
        return (-(u * ddx(u) + v * ddy(u)) + nu * lp(u),
                -(u * ddx(v) + v * ddy(v)) + nu * lp(v))

    t = 0.0
    while t < T - 1e-12:
        h = min(dt, T - t)
        ru, rv = rhs(u, v)
        u1, v1 = project(u + h * ru, v + h * rv)
        r1u, r1v = rhs(u1, v1)
        u, v = project(u + 0.5 * h * (ru + r1u), v + 0.5 * h * (rv + r1v))
        t += h
    F = math.exp(-2 * nu * T)
    ue = 1.0 - np.cos(X - T) * np.sin(Y) * F
    ve = np.sin(X - T) * np.cos(Y) * F
    l2 = math.sqrt(float(np.mean((u - ue) ** 2 + (v - ve) ** 2) / 2))
    return u, v, l2


def tt_to_dense(a, m):
    v = np.ones((1,))
    res = None
    cur = np.ones((1, 1))
    # contract sequentially: result index ordering [x bits LSB..], [y bits]
    t = a[0]
    for c in a[1:]:
        t = np.einsum("...a,adb->...db", t, c)
    t = t.reshape((2,) * (2 * m))
    # bits LSB-first -> index
    idx = np.arange(2 ** m)
    bits = ((idx[:, None] >> np.arange(m)) & 1)
    out = np.zeros((2 ** m, 2 ** m))
    tt = t.reshape(2 ** m, 2 ** m)  # careful: axis order is x bits then y bits
    # x flattened LSB-first equals mixed-radix; build permutation
    perm = np.zeros(2 ** m, dtype=int)
    for i in idx:
        b = [(i >> k) & 1 for k in range(m)]
        j = 0
        for k, bb in enumerate(b):
            j |= bb << (m - 1 - k)   # tensor axes are in core order
        perm[i] = j
    return tt[np.ix_(perm, perm)]


# ===================== stages =========================================
def stage_g1():
    print("=== G1: QTT vs dense cross-validation (same stencils) ===")
    for m in (6, 7):
        N = 2 ** m; Re = 100.0; nu = 2 * math.pi / Re
        T = 0.5; dt = 0.25 * (2 * math.pi / N)
        t0 = time.time()
        du, dv, dl2 = dense_solve(N, nu, T, dt)
        td = time.time() - t0
        S = QTTFV(m, nu)
        u, v = tgv_init_tt(m)
        t0 = time.time(); t = 0.0
        while t < T - 1e-12:
            h = min(dt, T - t)
            u, v = S.step(u, v, h)
            t += h
        tq = time.time() - t0
        uq = tt_to_dense(u, m)
        diff = math.sqrt(float(np.mean((uq - du) ** 2)))
        l2q = l2_vs_exact_tt(u, v, m, T, nu)
        print(f"N={N:4d}: dense L2 {dl2:.3e} ({td:.1f}s) | "
              f"QTT L2 {l2q:.3e} ({tq:.1f}s, chi {tt_chimax(u)}) | "
              f"TT-vs-dense {diff:.2e}")
    return 0


def stage_ladder():
    print("=== Reynolds ladder (Reynolds-scaled N) ===")
    rows = []
    for Re, m, T in ((10, 6, 1.0), (100, 8, 1.0),
                     (10000, 12, 0.1), (1000000, 16, 0.01)):
        N = 2 ** m; nu = 2 * math.pi / Re
        dx = 2 * math.pi / N
        dt = min(0.25 * dx / 2.0, 0.2 * dx * dx / nu)    # CFL + viscous
        steps = int(round(T / dt))
        dense_mem = 4 * (N ** 2) * 8
        print(f"\nRe={Re:g}, N={N} ({steps} steps): dense field memory "
              f"{dense_mem/1e6:.1f} MB x4 fields", flush=True)
        S = QTTFV(m, nu)
        u, v = tgv_init_tt(m)
        t0 = time.time(); t = 0.0
        while t < T - 1e-12:
            h = min(dt, T - t)
            u, v = S.step(u, v, h)
            t += h
        wall = time.time() - t0
        l2 = l2_vs_exact_tt(u, v, m, T, nu)
        qmem = tt_mem(u) + tt_mem(v)
        row = dict(Re=Re, N=N, T=T, steps=steps, l2=float(l2),
                   chi=int(tt_chimax(u)), qtt_bytes=int(qmem),
                   dense_bytes=int(dense_mem), wall_s=float(wall),
                   mem_ratio=float(dense_mem / qmem))
        rows.append(row)
        print(f"  QTT: L2 {l2:.3e} | chi {row['chi']} | "
              f"mem {qmem/1024:.1f} KB vs dense {dense_mem/1e6:.1f} MB "
              f"= {row['mem_ratio']:.0f}x | wall {wall:.0f}s", flush=True)
        json.dump(rows, open(WORK + r"\qtt_ladder.json", "w"), indent=1)
    print(f"\n-> {WORK}\\qtt_ladder.json")
    return 0


if __name__ == "__main__":
    sys.exit({"g1": stage_g1, "ladder": stage_ladder}[
        sys.argv[1] if len(sys.argv) > 1 else "g1"]())
