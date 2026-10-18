"""Dark-sector fifth force: DM feels G(1+2 beta^2) toward DM, G toward baryons; baryons feel gravity only.
LCDM background fixed at the Planck+DESI+SN best fit. Bounds on beta^2 from fsigma8 (Gold 2018)."""
import numpy as np
from scipy.integrate import solve_ivp
from growth import ZD, FS8, IC, _hda, HDA_REF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

H0, OM, OB = 68.2, 0.304, 0.049
OC = OM - OB
SIGMA8_PLANCK = 0.811            # LCDM amplitude at z = 0 implied by the CMB normalisation
E = lambda z: np.sqrt(OM * (1 + z) ** 3 + 1 - OM)


def growth(beta2, zi=50.0):
    """Returns (delta_c, delta_b, dlnD/dx for each) as functions of z, normalised to delta = a at z = zi."""
    def rhs(x, y):
        a = np.exp(x); z = 1 / a - 1
        E2 = E(z) ** 2
        dlnE = -1.5 * OM * a ** -3 / E2
        src = 1.5 / (E2 * a ** 3)
        dc, vc, db, vb = y
        return [vc, -(2 + dlnE) * vc + src * (OB * db + (1 + 2 * beta2) * OC * dc),
                vb, -(2 + dlnE) * vb + src * (OB * db + OC * dc)]
    ai = 1 / (1 + zi)
    s = solve_ivp(rhs, [np.log(ai), 0], [ai, ai, ai, ai], dense_output=True, rtol=1e-8, atol=1e-11)
    return s


def chi2(beta2, tracer="dm", normalise="marginalise"):
    s = growth(beta2)
    s0 = growth(0.0)
    i = 0 if tracer == "dm" else 2
    y0 = s.sol(0)
    nu = np.array([s.sol(np.log(1 / (1 + z)))[i + 1] for z in ZD]) / y0[i]       # f D / D(0), shape
    mu = nu / (_hda(E, 0.0, ZD) / HDA_REF)
    if normalise == "marginalise":
        s8 = (FS8 @ IC @ mu) / (mu @ IC @ mu)
    else:  # CMB-normalised amplitude: sigma8 today scales with the extra growth since z_i
        s8 = SIGMA8_PLANCK * y0[i] / s0.sol(0)[i]
    r = FS8 - s8 * mu
    return float(r @ IC @ r), float(s8)


betas = np.linspace(0, 0.4, 81)
lines = ["Dark-sector fifth force on linear growth (22 fsigma8 points, Gold 2018); LCDM background fixed", "",
         f"{'beta^2':>7s} {'DM tracer, s8 marg.':>22s} {'baryon tracer, s8 marg.':>25s} {'DM tracer, CMB-normalised':>27s} {'baryon, CMB-norm.':>19s}"]
res = {k: [] for k in ["dm_m", "b_m", "dm_c", "b_c"]}
for b2 in betas:
    c1, s1 = chi2(b2, "dm", "marginalise"); c2, s2 = chi2(b2, "baryon", "marginalise")
    c3, s3 = chi2(b2, "dm", "cmb"); c4, s4 = chi2(b2, "baryon", "cmb")
    for k, v in zip(res, [c1, c2, c3, c4]):
        res[k].append(v)
    if abs(b2 * 20 - round(b2 * 20)) < 1e-9:
        lines.append(f"{b2:7.3f} {c1:12.2f} (s8={s1:.3f}) {c2:15.2f} (s8={s2:.3f}) {c3:17.2f} (s8={s3:.3f}) {c4:9.2f} (s8={s4:.3f})")
lines.append("")
for k, lab in zip(res, ["DM tracer, sigma8 marginalised", "baryon tracer, sigma8 marginalised", "DM tracer, CMB-normalised sigma8", "baryon tracer, CMB-normalised sigma8"]):
    c = np.array(res[k]); c0 = c[0]
    above = betas[c - c0 > 3.84]
    bound = f"beta^2 < {above[0]:.3f} (95%)" if len(above) else "no 95% bound within beta^2 < 0.4"
    best = betas[np.argmin(c)]
    lines.append(f"  {lab:38s}: best beta^2 = {best:.3f}, chi2(0) = {c0:.2f}, min = {c.min():.2f};  {bound}")
text = "\n".join(lines); print(text)
open("results/growth_force.txt", "w").write(text + "\n")

fig, ax = plt.subplots(figsize=(7, 4.5))
for k, lab in zip(res, ["DM tracer, σ8 marginalised", "baryon tracer, σ8 marginalised", "DM tracer, CMB-normalised", "baryon tracer, CMB-normalised"]):
    ax.plot(betas, np.array(res[k]) - res[k][0], label=lab)
ax.axhline(3.84, color="gray", ls=":", label="Δχ² = 3.84 (95%)")
ax.set_xlabel("β² (fifth-force strength relative to gravity, DM–DM)"); ax.set_ylabel("Δχ² vs β = 0")
ax.set_ylim(-2, 30); ax.legend(fontsize=8)
fig.savefig("results/growth_force.png", dpi=120, bbox_inches="tight")
