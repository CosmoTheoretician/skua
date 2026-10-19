"""Aggregate halo baryon fractions over the production grid into f_b(beta)/f_b(0) on common mass bins."""
import glob, json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EDGES = np.array([1e-3, 2e-3, 4e-3, 8e-3, 1e-1])       # code mass; 1 code mass ~ 1.3e12 Msun/h for a 10 Mpc/h box
runs = {}
for cfgf in glob.glob(os.path.expanduser("~/scraper/sim_gpu2/out/*/config.json")):
    c = json.load(open(cfgf)); c = {**c.get("derived", {}), **c.get("config", c)}
    if not str(c.get("label", "")).startswith("prod_"):
        continue
    name = os.path.basename(os.path.dirname(cfgf)); f = f"results/halos_{name}.npy"
    if os.path.exists(f):
        runs[(round(float(c["fifth_lambda"]), 3), float(c["fifth_beta2"]))] = np.load(f)

def binned(rows, col):
    return np.array([np.nanmedian(rows[(rows[:, 0] >= lo) & (rows[:, 0] < hi), col]) for lo, hi in zip(EDGES[:-1], EDGES[1:])])

lines = ["Dark-sector fifth force: halo baryon fraction relative to the beta = 0 run (same seed), common mass bins", ""]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for lam in sorted({k[0] for k in runs}):
    base = runs.get((lam, 1e-9)) if (lam, 1e-9) in runs else runs.get((lam, 0.0))
    if base is None:
        continue
    b200, bin_ = binned(base, 2), binned(base, 3)
    lines += [f"range lambda = {lam} box units  (baseline f_b/f_cosmic at R200: " + ", ".join(f"{v:.3f}" for v in b200) + ")", "",
              f"{'beta^2':>7s} | " + " ".join(f"{'bin'+str(i+1)+' R200':>10s}" for i in range(4)) + " | " + " ".join(f"{'bin'+str(i+1)+' inner':>11s}" for i in range(4))]
    for (l, b2), rows in sorted(runs.items()):
        if l != lam or b2 < 1e-6:
            continue
        r200, rin = binned(rows, 2) / b200, binned(rows, 3) / bin_
        lines.append(f"{b2:7.2f} | " + " ".join(f"{v:10.3f}" for v in r200) + " | " + " ".join(f"{v:11.3f}" for v in rin))
        ax = axes[0 if lam == min({k[0] for k in runs}) else 1]
        ax.plot(range(1, 5), r200, "o-", label=f"β² = {b2}, R200"); ax.plot(range(1, 5), rin, "s--", label=f"β² = {b2}, inner")
    lines.append("")
lines += ["bins: M200 in code mass [1-2, 2-4, 4-8, 8-100] x 1e-3 (lightest to heaviest); inner = within R200/3.",
          "Halos here have R200 ~ 0.006-0.01 box units, i.e. lambda/R200 ~ 2 (0.016) and ~4 (0.03): the force spans the halo.",
          "Numerical baseline depletion (unequal particle masses) cancels in the ratio to first order; see the equal-mass check."]
text = "\n".join(lines); print(text)
open("results/fifth_force_map.txt", "w").write(text + "\n")
for ax, lam in zip(axes, sorted({k[0] for k in runs})):
    ax.axhline(1, color="gray", lw=0.8); ax.set_title(f"range λ = {lam} box units"); ax.set_xlabel("mass bin (light → heavy)")
    ax.set_xticks(range(1, 5)); ax.legend(fontsize=7)
axes[0].set_ylabel("f_b(β) / f_b(0)")
fig.savefig("results/fifth_force_map.png", dpi=120, bbox_inches="tight")
