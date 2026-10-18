"""Finite-range dark fifth force: G_eff(k) = G [1 + 2 beta^2 k^2/(k^2 + a^2/lambda^2)] for DM-DM.
Bound on beta^2 from fsigma8 (CMB-normalised amplitude) as a function of the range lambda."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from growth import ZD, FS8, IC, _hda, HDA_REF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OM, OB = 0.304, 0.049
OC = OM - OB
SIGMA8_PLANCK = 0.811
E = lambda z: np.sqrt(OM * (1 + z) ** 3 + 1 - OM)


def growth(beta2, lam, k, zi=50.0):
    def rhs(x, y):
        a = np.exp(x); z = 1 / a - 1
        E2 = E(z) ** 2
        dlnE = -1.5 * OM * a ** -3 / E2
        src = 1.5 / (E2 * a ** 3)
        yuk = k**2 / (k**2 + (a / lam) ** 2)
        dc, vc, db, vb = y
        return [vc, -(2 + dlnE) * vc + src * (OB * db + (1 + 2 * beta2 * yuk) * OC * dc),
                vb, -(2 + dlnE) * vb + src * (OB * db + OC * dc)]
    ai = 1 / (1 + zi)
    return solve_ivp(rhs, [np.log(ai), 0], [ai, ai, ai, ai], dense_output=True, rtol=1e-8, atol=1e-11)


def chi2(beta2, lam, k, tracer_index=0):
    s, s0 = growth(beta2, lam, k), growth(0.0, lam, k)
    y0 = s.sol(0)
    nu = np.array([s.sol(np.log(1 / (1 + z)))[tracer_index + 1] for z in ZD]) / y0[tracer_index]
    mu = nu / (_hda(E, 0.0, ZD) / HDA_REF)
    s8 = SIGMA8_PLANCK * y0[tracer_index] / s0.sol(0)[tracer_index]
    r = FS8 - s8 * mu
    return float(r @ IC @ r)


def bound(lam, k, cap=3.0):
    c0 = chi2(0.0, lam, k)
    f = lambda b2: chi2(b2, lam, k) - c0 - 3.84
    return brentq(f, 0.0, cap) if f(cap) > 0 else np.inf


lams = np.logspace(-1, 3, 25)            # Mpc/h
ks = [0.05, 0.15, 0.3]                   # h/Mpc: the range of scales fsigma8 measurements probe
lines = ["95% upper bound on beta^2 (DM-DM fifth force) vs range lambda, fsigma8 with CMB-normalised sigma8", "",
         f"{'lambda [Mpc/h]':>15s}" + "".join(f" {'k='+str(k)+' h/Mpc':>16s}" for k in ks)]
B = {k: [] for k in ks}
for lam in lams:
    row = f"{lam:15.2f}"
    for k in ks:
        b = bound(lam, k); B[k].append(b)
        row += f" {('> 3 (unbounded)' if not np.isfinite(b) else f'{b:.4f}'):>16s}"
    lines.append(row)
lines += ["", "Reading: the linear bound is the long-range value (~0.01) for lambda >> 1/k, weakens as lambda drops below",
          "the probed scales, and disappears once the force is shorter than a few Mpc/h. Below that, only nonlinear",
          "structure (halos, clusters, galaxies) can constrain it."]
text = "\n".join(lines); print(text)
open("results/growth_force_range.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(7, 4.5))
for k in ks:
    b = np.array(B[k]); ax.plot(lams, np.where(np.isfinite(b), b, 3.0), label=f"k = {k} h/Mpc")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_ylim(5e-3, 3)
ax.set_xlabel("force range λ [Mpc/h]"); ax.set_ylabel("95% bound on β²"); ax.legend()
ax.axhline(1, color="gray", ls=":"); ax.text(0.12, 1.1, "force as strong as gravity", fontsize=8, color="gray")
fig.savefig("results/growth_force_range.png", dpi=120, bbox_inches="tight")
