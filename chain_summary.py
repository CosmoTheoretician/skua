"""Means, 68% intervals and (w0, wa) contours from the cobaya chains that exist so far."""
import glob, os
import numpy as np
from getdist import loadMCSamples, plots
import matplotlib
matplotlib.use("Agg")

runs = [("mcmc_w0wa", "Planck + DESI DR2 + DES-SN5YR (2024)"), ("mcmc_w0wa_pp", "Planck + DESI DR2 + Pantheon+"),
        ("mcmc_w0wa_dovekie", "Planck + DESI DR2 + DES-Dovekie")]
samples, labels = [], []
for name, lab in runs:
    if not glob.glob(f"chains/{name}.1.txt"):
        continue
    s = loadMCSamples(f"chains/{name}", settings={"ignore_rows": 0.3})
    samples.append(s); labels.append(lab)
    n = s.numrows
    print(f"\n{lab}  ({n} samples after burn-in)")
    for p in ["w", "wa", "H0", "omegam", "sigma8", "rdrag"]:
        try:
            m, (lo, hi) = s.mean(p), s.twoTailLimits(p, 0.68)
            print(f"  {p:7s} = {m:8.4f}  (68%: {lo:.4f} to {hi:.4f})")
        except Exception:
            pass
    try:
        print("  R-1 =", round(s.getGelmanRubin(), 4))
    except Exception:
        pass
if samples:
    g = plots.get_subplot_plotter(width_inch=5)
    g.plot_2d(samples, "w", "wa", filled=True, lims=[-1.3, -0.4, -2.5, 0.5])
    g.add_x_marker(-1); g.add_y_marker(0); g.add_legend(labels, fontsize=8)
    g.export("results/w0wa_contours.png")
    print("\nfigure: results/w0wa_contours.png")
