import os, sys, time, json, math
os.environ.setdefault("OMP_NUM_THREADS", "4")
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qtt_fv as Q
Re, m, T = 100, 8, 1.0
N = 2 ** m; nu = 2 * math.pi / Re; dx = 2 * math.pi / N
dt = min(0.25 * dx / 2.0, 0.2 * dx * dx / nu); steps = int(round(T / dt))
S = Q.QTTFV(m, nu); u, v = Q.tgv_init_tt(m); t = 0.0; chis = []; t0 = time.time()
while t < T - 1e-12:
    h = min(dt, T - t); u, v = S.step(u, v, h); t += h; chis.append(int(Q.tt_chimax(u)))
    if len(chis) % 100 == 0: print(f"  step {len(chis)}/{steps}  chi {chis[-1]}  {time.time()-t0:.0f}s", flush=True)
wall = time.time() - t0; l2 = float(Q.l2_vs_exact_tt(u, v, m, T, nu))
rec = dict(Re=Re, N=N, T=T, dt=dt, steps=steps, l2=l2, chi_max=max(chis), qtt_bytes=int(Q.tt_mem(u) + Q.tt_mem(v)),
           dense_bytes=int(4 * N * N * 8), wall_s=round(wall, 1), solver="qtt_fv.py with gesvd fallback, 2026-09-14",
           ladder_cited=dict(l2=6.179507392404894e-05, qtt_bytes=13472, wall_s=3515.1))
json.dump(rec, open(os.path.join(HERE, "rerun_Re100_20260914.json"), "w"), indent=1); print("done", rec, flush=True)
