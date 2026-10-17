"""fsigma8 likelihood (Gold 2018 compilation, 22 points) for an arbitrary background."""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

D = np.array([
    [0.02, 0.428, 0.0465, 0.3], [0.02, 0.398, 0.065, 0.3], [0.02, 0.314, 0.048, 0.266], [0.10, 0.370, 0.130, 0.3],
    [0.15, 0.490, 0.145, 0.31], [0.17, 0.510, 0.060, 0.3], [0.18, 0.360, 0.090, 0.27], [0.38, 0.440, 0.060, 0.27],
    [0.25, 0.3512, 0.0583, 0.25], [0.37, 0.4602, 0.0378, 0.25], [0.32, 0.384, 0.095, 0.274], [0.59, 0.488, 0.060, 0.307115],
    [0.44, 0.413, 0.080, 0.27], [0.60, 0.390, 0.063, 0.27], [0.73, 0.437, 0.072, 0.27], [0.60, 0.550, 0.120, 0.3],
    [0.86, 0.400, 0.110, 0.3], [1.40, 0.482, 0.116, 0.27], [0.978, 0.379, 0.176, 0.31], [1.23, 0.385, 0.099, 0.31],
    [1.526, 0.342, 0.070, 0.31], [1.944, 0.364, 0.106, 0.31]])
ZD, FS8, SIG, OMREF = D.T
C = np.diag(SIG**2)
C[12:15, 12:15] = 1e-3 * np.array([[6.400, 2.570, 0.000], [2.570, 3.969, 2.540], [0.000, 2.540, 5.184]])
C[18:22, 18:22] = 1e-2 * np.array([[3.098, 0.892, 0.329, -0.021], [0.892, 0.980, 0.436, 0.076],
                                   [0.329, 0.436, 0.490, 0.350], [-0.021, 0.076, 0.350, 1.124]])
IC = np.linalg.inv(C)
N_GROWTH = len(ZD)


def _hda(E, ok, z):
    zg = np.linspace(0, z.max() * 1.001, 2001)
    chi = cumulative_trapezoid(1 / E(zg), zg, initial=0)
    dm = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok) if ok > 1e-8 else chi
    return E(z) * np.interp(z, zg, dm) / (1 + z)


_E_ref = lambda om: (lambda z: np.sqrt(om * (1 + z) ** 3 + 1 - om))
HDA_REF = np.array([_hda(_E_ref(o), 0.0, np.array([z]))[0] for z, o in zip(ZD, OMREF)])


def chi2_growth(E, rhom, ok, zi=50.0):
    """E(z): expansion rate / H0.  rhom(z): matter density in units of rho_crit,0.  ok: Omega_k.
    GR growth equation with sigma8 marginalised analytically. Returns (chi2, sigma8_best)."""
    def dlnE(x, h=1e-4):
        z1, z2 = np.exp(-(x + h)) - 1, np.exp(-(x - h)) - 1
        return (np.log(E(z1)) - np.log(E(z2))) / (2 * h)

    def rhs(x, y):
        z = np.exp(-x) - 1
        return [y[1], -(2 + dlnE(x)) * y[1] + 1.5 * rhom(z) / E(z) ** 2 * y[0]]

    ai = 1 / (1 + zi)
    s = solve_ivp(rhs, [np.log(ai), 0], [ai, ai], dense_output=True, rtol=1e-7, atol=1e-10)
    d0 = s.sol(0)[0]
    nu = np.array([s.sol(np.log(1 / (1 + z)))[1] / d0 for z in ZD])        # f*D/D0
    mu = nu / (_hda(E, ok, ZD) / HDA_REF)                                   # AP-corrected model shape
    s8 = (FS8 @ IC @ mu) / (mu @ IC @ mu)
    r = FS8 - s8 * mu
    return float(r @ IC @ r), float(s8)
