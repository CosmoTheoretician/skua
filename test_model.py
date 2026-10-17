"""python test_model.py <model> [pantheon|des] [cmb]
A model file defines NAMES, BOUNDS, P0, E(z,p), OK(p) and optionally OM(p) (total matter today; default p[0])."""
import sys, os, json, importlib
import numpy as np
from scipy.optimize import differential_evolution, minimize
from likelihoods import SN, chi2_cc, chi2_bao, ZG
from cmb import chi2_cmb

name = sys.argv[1]
which = sys.argv[2] if len(sys.argv) > 2 else "pantheon"
use_cmb = "cmb" in sys.argv[3:]
tag = which + ("_cmb" if use_cmb else "")
m = importlib.import_module(f"models.{name}")
OM = getattr(m, "OM", lambda p: p[0])
sn = SN(which)
REFFILE = f"results/ref_{tag}.json"
REF = json.load(open(REFFILE)) if os.path.exists(REFFILE) else {}


def parts(q):
    if use_cmb:
        H0, wb, p = q[0], q[1], q[2:]
    else:
        H0, rd, p = q[0], q[1], q[2:]
    E = lambda z: m.E(z, p)
    with np.errstate(all="ignore"):
        e = E(ZG)
    if not np.all(np.isfinite(e)) or np.any(e <= 0):
        return None
    ok = OM and m.OK(p)
    out = {}
    if use_cmb:
        c, rd, R, lA = chi2_cmb(E, ok, OM(p), H0 / 100, wb)
        if not np.isfinite(c):
            return None
        out["CMB"] = c; out["_rd"] = rd; out["_R"] = R; out["_lA"] = lA
    out["CC"] = chi2_cc(E, H0); out["BAO"] = chi2_bao(E, ok, H0, rd); out["SN"] = sn.chi2(E, ok)
    return out


def total(q):
    o = parts(q)
    return 1e9 if o is None else sum(v for k, v in o.items() if not k.startswith("_"))


bounds = [(50, 90), (0.019, 0.026) if use_cmb else (100, 200)] + list(m.BOUNDS)
de = differential_evolution(total, bounds, seed=1, popsize=20, maxiter=300, tol=1e-8, polish=False)
res = minimize(total, de.x, method="Nelder-Mead", options=dict(maxiter=20000, xatol=1e-7, fatol=1e-7))
q, c = (res.x, res.fun) if res.fun < de.fun else (de.x, de.fun)
o = parts(q); k = len(q)
pnames = ["H0", "wb" if use_cmb else "rd"] + list(m.NAMES)
lines = [f"model: {name}   data: CC + BAO + {sn.name}" + (" + Planck distance priors" if use_cmb else "") + f"   parameters: " + ", ".join(pnames), "",
         f"chi2 = {c:.2f}  (" + ", ".join(f"{kk} {v:.2f}" for kk, v in o.items() if not kk.startswith("_")) + f");  k = {k};  AIC = {c + 2*k:.2f}"]
if "flat LCDM" in REF:
    lines.append(f"vs flat LCDM: dchi2 = {c - REF['flat LCDM']:+.1f}, dAIC = {c + 2*k - REF['flat LCDM'] - 2*REF['k_lcdm']:+.1f}")
lines += ["", "best fit: " + ", ".join(f"{n} = {v:.4g}" for n, v in zip(pnames, q))]
if use_cmb:
    lines.append(f"derived: r_d = {o['_rd']:.1f} Mpc, R = {o['_R']:.4f}, l_A = {o['_lA']:.2f}")
if name == "lcdm":
    REF.update({"flat LCDM": c, "k_lcdm": k}); json.dump(REF, open(REFFILE, "w"))
text = "\n".join(lines); print(text)
open(f"results/{name}_{tag}.txt", "w").write(text + "\n")
