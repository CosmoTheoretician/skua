"""Box-mode growth in a two-species run vs linear two-fluid theory, for either the Janus coupling or the
fifth-force mode (A18), static or expanding. No free parameters: seeds are Poisson shot noise, 1/N per mode."""
import sys, json, glob, os
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RUN = sys.argv[1]
cfg = json.load(open(RUN + "/config.json")); cfg = {**cfg.get("derived", {}), **cfg.get("config", cfg)}
files = sorted(glob.glob(RUN + "/pts_*.npz")); d0 = np.load(files[0])
box = float(d0["box"]); N = int(cfg["npos"]); ratio = float(cfg["mass_ratio"])
b2 = float(cfg.get("fifth_beta2", 0.0)); lam = float(cfg.get("fifth_lambda", 0.0)) * box
fifth = b2 > 0
TJ = float(cfg["t_ref"]); expansion = bool(cfg.get("expansion", False)); T0 = float(cfg.get("t0_exp", 0.1)) * TJ
G = float(cfg["G"])
G_RHO_NEG = G * int(cfg["nneg"]) * float(cfg["w_neg"]) / box**3    # from the run's own G and weights
G_RHO_POS = G * N * float(cfg["w_pos"]) / box**3
print(f"G rho_neg = {G_RHO_NEG:.3f}, G rho_pos = {G_RHO_POS:.3f}, t_J = {TJ:.4f}")
k1 = 2 * np.pi / box
yuk = lambda k: 1.0 if lam <= 0 else k**2 / (k**2 + 1 / lam**2)
print(f"run {os.path.basename(RUN)}: mode = {'fifth force' if fifth else 'Janus'}, ratio = {ratio}, beta2 = {b2}, "
      f"lambda = {lam:.4g}, expansion = {expansion}, N = {N}")

KV = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, -1, 0], [1, 0, 1], [1, 0, -1], [0, 1, 1], [0, 1, -1]], float)


def power(x):
    ph = 2 * np.pi * (x / box) @ KV.T
    return np.abs(np.exp(-1j * ph).mean(axis=0)) ** 2 - 1.0 / len(x)


ts, Pp, Pn = [], [], []
for f in files:
    d = np.load(f); ts.append(float(d["t"]))
    Pp.append(power(d["pos"].astype(float)).mean()); Pn.append(power(d["neg"].astype(float)).mean())
ts, Pp, Pn = map(np.array, (ts, Pp, Pn))

a_of = (lambda t: 1.0 if (not expansion or t < T0) else (t / T0) ** (2 / 3))
H_of = (lambda t: 0.0 if (not expansion or t < T0) else 2 / (3 * t))
# coupling matrix M: d2 delta_i = (4 pi / a^3) sum_j M_ij G rho_j delta_j ; i,j in (pos, neg)
if fifth:   # both positive, neg = DM with Yukawa self-attraction at the box scale (|k| = k1 and sqrt2 k1 averaged)
    yb = 0.5 * (yuk(k1) + yuk(np.sqrt(2) * k1))
    M = np.array([[1.0, 1.0], [1.0, 1.0 + 2 * b2 * yb]])
else:       # Janus: like attracts, unlike repels
    M = np.array([[1.0, -1.0], [-1.0, 1.0]])
GR = np.array([G_RHO_POS, G_RHO_NEG])


def rhs(t, y):
    a, H = a_of(t), H_of(t)
    dp, vp, dn, vn = y
    src = 4 * np.pi / a**3
    return [vp, -2 * H * vp + src * (M[0, 0] * GR[0] * dp + M[0, 1] * GR[1] * dn),
            vn, -2 * H * vn + src * (M[1, 0] * GR[0] * dp + M[1, 1] * GR[1] * dn)]


tt = np.linspace(1e-6, ts.max(), 1500)
resp = {n: solve_ivp(rhs, [tt[0], tt[-1]], y0, t_eval=tt, rtol=1e-9, atol=1e-12, max_step=1e-3).y
        for n, y0 in [("p", [1, 0, 0, 0]), ("n", [0, 0, 1, 0])]}
P0 = 1.0 / N
Pp_th = P0 * (resp["p"][0] ** 2 + resp["n"][0] ** 2)
Pn_th = P0 * (resp["p"][2] ** 2 + resp["n"][2] ** 2)
noise = 1 / np.sqrt(len(d0["pos"]))

lines = [f"{'t/tJ':>6s} {'neg(DM): sim':>13s} {'theory':>9s} {'ratio':>6s}   {'pos(baryon): sim':>16s} {'theory':>9s} {'ratio':>6s}"]
for t, pp, pn in zip(ts, Pp, Pn):
    tn, tp = np.sqrt(np.interp(t, tt, Pn_th)), np.sqrt(np.interp(t, tt, Pp_th))
    sn, sp = np.sqrt(max(pn, 0)), np.sqrt(max(pp, 0))
    flag = "  (nonlinear)" if tn > 0.3 else ("  (noise)" if sn < 1.5 * noise else "")
    lines.append(f"{t/TJ:6.2f} {sn:13.2e} {tn:9.2e} {sn/tn:6.2f}   {sp:16.2e} {tp:9.2e} {sp/tp:6.2f}{flag}")
text = "\n".join(lines); print(text)
open(f"results/nbody_growth_{os.path.basename(RUN)}.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(tt / TJ, np.sqrt(Pn_th), "C1", label="DM / neg, linear theory"); ax.plot(tt / TJ, np.sqrt(Pp_th), "C0", label="baryons / pos, linear theory")
ax.plot(ts / TJ, np.sqrt(np.clip(Pn, 0, None)), "C1o", ms=4, label="N-body"); ax.plot(ts / TJ, np.sqrt(np.clip(Pp, 0, None)), "C0s", ms=4)
ax.axhline(noise, color="gray", ls=":", label="subsample noise"); ax.axhline(1, color="k", lw=0.5)
ax.set_yscale("log"); ax.set_xlabel("t / t_J"); ax.set_ylabel("δ_rms, lowest box modes"); ax.legend(fontsize=8)
ax.set_title(f"{os.path.basename(RUN)}: {'fifth force β²=' + str(b2) if fifth else 'Janus'}, ratio {ratio}")
fig.savefig(f"results/nbody_growth_{os.path.basename(RUN)}.png", dpi=120, bbox_inches="tight")
