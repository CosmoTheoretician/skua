"""Is the evolving-DE signal a low-redshift supernova calibration offset? Profile a free magnitude
offset Delta for z < z_split alongside the absolute magnitude, and see what survives."""
import sys, importlib
import numpy as np
from scipy.optimize import differential_evolution, minimize
from likelihoods import SN, chi2_cc, chi2_bao, ZG, comoving, transverse
from cmb import chi2_cmb

Z_SPLIT = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1


class SNoff(SN):
    def __init__(self, which):
        super().__init__(which)
        self.A = np.column_stack([self.ones, (self.z < Z_SPLIT).astype(float)])
        self.AtC = self.A.T @ self.icov
        self.M = np.linalg.inv(self.AtC @ self.A)
        self.nlow = int((self.z < Z_SPLIT).sum())

    def chi2_off(self, E, ok):
        chi = comoving(E)
        if chi is None:
            return np.inf, 0.0
        with np.errstate(all="ignore"):
            mu = 5 * np.log10((1 + self.zhel) * np.interp(self.z, ZG, transverse(chi, ok)))
        if not np.all(np.isfinite(mu)):
            return np.inf, 0.0
        r = self.mu - mu
        b = self.AtC @ r
        x = self.M @ b                      # (M, Delta) at the minimum
        return float(r @ self.icov @ r - b @ x), float(x[1])


def fit(model, sn, use_off):
    m = importlib.import_module(f"models.{model}")
    OM = getattr(m, "OM", lambda p: p[0])

    def total(q):
        H0, wb, p = q[0], q[1], q[2:]
        E = lambda z: m.E(z, p)
        with np.errstate(all="ignore"):
            e = E(ZG)
        if not np.all(np.isfinite(e)) or np.any(e <= 0):
            return 1e9
        ok = m.OK(p)
        c_cmb = chi2_cmb(E, ok, OM(p), H0 / 100, wb)[0]
        c_sn = sn.chi2_off(E, ok)[0] if use_off else sn.chi2(E, ok)
        return c_cmb + chi2_cc(E, H0) + chi2_bao(E, ok, H0, 147.0 * 0 + rd_of(E, ok, OM(p), H0 / 100, wb)) + c_sn

    def rd_of(E, ok, Om, h, wb):
        return chi2_cmb(E, ok, Om, h, wb)[1]

    bounds = [(55, 80), (0.020, 0.025)] + list(m.BOUNDS)
    de = differential_evolution(total, bounds, seed=1, popsize=15, maxiter=200, tol=1e-7, polish=False)
    res = minimize(total, de.x, method="Nelder-Mead", options=dict(maxiter=10000, xatol=1e-6, fatol=1e-6))
    q, c = (res.x, res.fun) if res.fun < de.fun else (de.x, de.fun)
    p = q[2:]; E = lambda z: m.E(z, p)
    delta = sn.chi2_off(E, m.OK(p))[1] if use_off else 0.0
    return c, q, delta


lines = [f"Low-z calibration offset test: free magnitude offset Delta for z < {Z_SPLIT} (profiled with M)", ""]
for which in ["des", "pantheon"]:
    sn = SNoff(which)
    lines.append(f"{sn.name}: {sn.nlow} SNe below z_split")
    res = {}
    for model in ["lcdm", "cpl"]:
        for off in [False, True]:
            c, q, d = fit(model, sn, off)
            res[(model, off)] = (c, q, d)
            tag = f"{model:5s} {'with Delta   ' if off else 'no offset    '}"
            pars = f"H0 = {q[0]:.1f}, Om = {q[2]:.3f}" + (f", w0 = {q[3]:.3f}, wa = {q[4]:.3f}" if model == "cpl" else "")
            lines.append(f"  {tag} chi2 = {c:9.2f}   {pars}" + (f"   Delta = {d:+.4f} mag" if off else ""))
    d0 = res[("cpl", False)][0] - res[("lcdm", False)][0]
    d1 = res[("cpl", True)][0] - res[("lcdm", True)][0]
    gain = res[("lcdm", False)][0] - res[("lcdm", True)][0]
    lines += [f"  CPL - LCDM:  {d0:+.1f} without offset  ->  {d1:+.1f} with offset   (signal surviving: {100 * d1 / d0:.0f}%)",
              f"  LCDM gains {gain:.1f} from the offset alone; offset LCDM wants: {res[('lcdm', True)][2]:+.4f} mag", ""]
lines += ["Context: survey-to-survey calibration systematics in current compilations are at the ~0.01-0.02 mag level;",
          "Efstathiou (2024) pointed to a ~0.04 mag offset between Pantheon+ and DES-SN5YR at z < 0.1."]
text = "\n".join(lines); print(text)
open(f"results/calib_test_z{Z_SPLIT}.txt", "w").write(text + "\n")
