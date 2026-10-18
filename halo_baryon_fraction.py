"""Halo baryon fractions from a two-species run: FoF on dark matter, R200 on total mass, f_b/f_cosmic per mass bin."""
import sys, os, json
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

RUN = sys.argv[1]; MIN_PART = int(sys.argv[2]) if len(sys.argv) > 2 else 200
d = np.load(RUN + "/final_state.npz", allow_pickle=True)
print("final_state keys:", {k: d[k].shape for k in d.files})
cfg = json.load(open(RUN + "/config.json")); cfg = {**cfg.get("derived", {}), **cfg.get("config", cfg)}
box = float(cfg["box"])
pos = d["pos"] if "pos" in d.files else d[[k for k in d.files if d[k].ndim == 2][0]]
sp = d["species"] if "species" in d.files else (d["s"] if "s" in d.files else None)
w = d["w"] if "w" in d.files else None
if sp is None:
    sys.exit("no species/sign array in final_state.npz; see keys above")
dm = sp < 0 if sp.min() < 0 else sp == 1            # 'neg' species = dark matter in fifth-force mode
xb, xd = pos[~dm] % box, pos[dm] % box
wb, wd = (w[~dm], w[dm]) if w is not None else (np.full((~dm).sum(), float(cfg["w_pos"])), np.full(dm.sum(), float(cfg["w_neg"])))
f_cosmic = wb.sum() / (wb.sum() + wd.sum())
print(f"baryons {len(xb)}, DM {len(xd)}, cosmic baryon fraction {f_cosmic:.4f}")

spacing = box / len(xd) ** (1 / 3); link = 0.2 * spacing
tree_d = cKDTree(xd, boxsize=box)
pairs = tree_d.query_pairs(link, output_type="ndarray")
n = len(xd)
g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
ncomp, lab = connected_components(g, directed=False)
counts = np.bincount(lab)
halos = np.where(counts >= MIN_PART)[0]
print(f"FoF: {ncomp} groups, {len(halos)} with >= {MIN_PART} DM particles")

tree_b = cKDTree(xb, boxsize=box)
rho_mean = (wb.sum() + wd.sum()) / box**3
rows = []
for h in halos:
    x = xd[lab == h]
    c = np.array([np.angle(np.exp(2j * np.pi * x[:, i] / box).mean()) for i in range(3)]) * box / (2 * np.pi) % box
    r200 = None
    for r in np.linspace(0.2 * link, 50 * link, 200)[::-1]:
        md = wd[tree_d.query_ball_point(c, r)].sum(); mb = wb[tree_b.query_ball_point(c, r)].sum()
        if (md + mb) / (4 / 3 * np.pi * r**3) >= 200 * rho_mean:
            r200, m200, mb200 = r, md + mb, mb
            break
    if r200 is not None:
        ri = r200 / 3
        mdi = wd[tree_d.query_ball_point(c, ri)].sum(); mbi = wb[tree_b.query_ball_point(c, ri)].sum()
        rows.append((m200, r200, mb200 / m200 / f_cosmic, mbi / (mdi + mbi) / f_cosmic))
rows = np.array(rows)
print(f"{len(rows)} halos with an R200")
edges = np.quantile(rows[:, 0], [0, 0.25, 0.5, 0.75, 1.0])
lines = [f"run {os.path.basename(RUN)}: beta2 = {cfg.get('fifth_beta2', 0)}, lambda = {cfg.get('fifth_lambda', 0)} box units", "",
         f"{'M200 bin (code mass)':>30s} {'n':>5s} {'f_b/f_cosmic (R200)':>20s} {'16-84%':>14s} {'inner (R200/3)':>15s}"]
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (rows[:, 0] >= lo) & (rows[:, 0] <= hi)
    q = np.percentile(rows[m, 2], [16, 50, 84]); qi = np.nanmedian(rows[m, 3])
    lines.append(f"{lo:14.3e} - {hi:12.3e} {m.sum():5d} {q[1]:20.3f} {q[0]:6.3f}-{q[2]:6.3f} {qi:15.3f}")
lines.append(f"{'all':>30s} {len(rows):5d} {np.median(rows[:, 2]):20.3f} {'':14s} {np.nanmedian(rows[:, 3]):15.3f}")
text = "\n".join(lines); print(text)
np.save(f"results/halos_{os.path.basename(RUN)}.npy", rows)
open(f"results/halos_{os.path.basename(RUN)}.txt", "w").write(text + "\n")
