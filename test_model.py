"""Fit any background model in models/<name>.py to CC + DESI BAO + Pantheon+ with r_d free.
A model file defines: NAMES (list), BOUNDS (list of (lo,hi)), P0 (list), E(z, p) -> E(z), OK(p) -> Omega_k."""
import sys, importlib
import numpy as np
from scipy.optimize import differential_evolution, minimize
from likelihoods import SN, chi2_cc, chi2_bao, ZG

REF = {"flat LCDM": 1422.50, "Janus": 1739.05}
name = sys.argv[1]
m = importlib.import_module(f"models.{name}")
sn = SN(sys.argv[2] if len(sys.argv) > 2 else "pantheon")


def total(q):
    H0, rd, p = q[0], q[1], q[2:]
    E = lambda z: m.E(z, p)
    with np.errstate(all="ignore"):
        e = E(ZG)
    if not np.all(np.isfinite(e)) or np.any(e <= 0):
        return 1e9
    ok = m.OK(p)
    return chi2_cc(E, H0) + chi2_bao(E, ok, H0, rd) + sn.chi2(E, ok)


bounds = [(50, 90), (100, 200)] + list(m.BOUNDS)
de = differential_evolution(total, bounds, seed=1, popsize=20, maxiter=300, tol=1e-8, polish=False)
res = minimize(total, de.x, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-7))
q, c = (res.x, res.fun) if res.fun < de.fun else (de.x, de.fun)
k = len(q)
E = lambda z: m.E(z, q[2:])
parts = chi2_cc(E, q[0]), chi2_bao(E, m.OK(q[2:]), q[0], q[1]), sn.chi2(E, m.OK(q[2:]))
lines = [f"model: {name}   parameters: H0, rd, " + ", ".join(m.NAMES), "",
         f"chi2 = {c:.2f}  (CC {parts[0]:.2f}, BAO {parts[1]:.2f}, SN {parts[2]:.2f});  k = {k};  AIC = {c + 2*k:.2f}",
         f"vs flat LCDM (k=3, AIC {REF['flat LCDM'] + 6:.2f}): dchi2 = {c - REF['flat LCDM']:+.1f}, dAIC = {c + 2*k - REF['flat LCDM'] - 6:+.1f}",
         f"vs Janus: dchi2 = {c - REF['Janus']:+.1f}", "",
         "best fit: " + ", ".join(f"{n} = {v:.4g}" for n, v in zip(["H0", "rd"] + list(m.NAMES), q))]
text = "\n".join(lines)
print(text)
open(f"results/{name}_{sys.argv[2] if len(sys.argv) > 2 else 'pantheon'}.txt", "w").write(text + "\n")
