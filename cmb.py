"""Compressed Planck 2018 CMB likelihood (R, l_A, omega_b) and the sound horizon from omega_b, omega_m."""
import numpy as np
from scipy.integrate import quad

C = 299792.458
OMEGA_G = 2.4728e-5                       # photons, T = 2.7255 K, in units of h^-2
OMEGA_R = OMEGA_G * (1 + 0.2271 * 3.046)  # + massless neutrinos
# Planck 2018 TT,TE,EE+lowE distance priors (Chen, Huang & Wang 2019): R, l_A, omega_b
MEAN = np.array([1.7502, 301.471, 0.02236])
SIG = np.array([0.0046, 0.090, 0.00015])
CORR = np.array([[1.0, 0.46, -0.66], [0.46, 1.0, -0.33], [-0.66, -0.33, 1.0]])
ICOV = np.linalg.inv(CORR * np.outer(SIG, SIG))


def z_star(wb, wm):  # Hu & Sugiyama 1996
    g1 = 0.0783 * wb ** -0.238 / (1 + 39.5 * wb ** 0.763)
    g2 = 0.560 / (1 + 21.1 * wb ** 1.81)
    return 1048 * (1 + 0.00124 * wb ** -0.738) * (1 + g1 * wm ** g2)


def z_drag(wb, wm):  # Eisenstein & Hu 1998
    b1 = 0.313 * wm ** -0.419 * (1 + 0.607 * wm ** 0.674)
    b2 = 0.238 * wm ** 0.223
    return 1291 * wm ** 0.251 / (1 + 0.659 * wm ** 0.828) * (1 + b1 * wb ** b2)


def with_radiation(E, h):
    orad = OMEGA_R / h**2
    return lambda z: np.sqrt(E(z) ** 2 + orad * ((1 + z) ** 4 - 1))


def sound_horizon(Ef, h, wb, z):
    """r_s(z) in Mpc: integral of c_s / H from z to infinity, in x = ln(1+z)."""
    def f(x):
        zz = np.exp(x) - 1
        cs = 1 / np.sqrt(3 * (1 + 3 * wb / (4 * OMEGA_G) / (1 + zz)))
        return cs / Ef(zz) * np.exp(x)
    return C / (100 * h) * quad(f, np.log(1 + z), np.log(1e9), limit=400)[0]


def comoving_transverse(Ef, h, ok, z):
    chi = quad(lambda zz: 1 / Ef(zz), 0, z, limit=400)[0]
    if ok > 1e-8:
        chi = np.sinh(np.sqrt(ok) * chi) / np.sqrt(ok)
    elif ok < -1e-8:
        chi = np.sin(np.sqrt(-ok) * chi) / np.sqrt(-ok)
    return C / (100 * h) * chi


def chi2_cmb_integral(E, ok, Om, h, wb):
    """older version: r_d from the integral (biased ~+1.7%); kept for reference."""
    wm = Om * h**2
    Ef = with_radiation(E, h)
    zs = z_star(wb, wm)
    dm = comoving_transverse(Ef, h, ok, zs)
    rs = sound_horizon(Ef, h, wb, zs)
    R = np.sqrt(Om) * (100 * h) * dm / C
    lA = np.pi * dm / rs
    d = np.array([R, lA, wb]) - MEAN
    rd = sound_horizon(Ef, h, wb, z_drag(wb, wm))
    return float(d @ ICOV @ d), rd, R, lA


def rd_aubourg(wb, wm, wnu=0.00064):
    """Aubourg et al. 2015 fit to CAMB, accurate to ~0.1%; wnu = 0.00064 is one 0.06 eV neutrino."""
    wcb = wm - wnu
    return 55.154 * np.exp(-72.3 * (wnu + 0.0006) ** 2) / (wcb ** 0.25351 * wb ** 0.12807)


def chi2_cmb(E, ok, Om, h, wb):
    """E: model E(z) without radiation. Returns (chi2, r_d, R, l_A). r_d from the Aubourg fit;
    r_s(z*) = r_d * (integral ratio), so the integral's ~1.7% bias cancels."""
    wm = Om * h**2
    Ef = with_radiation(E, h)
    zs = z_star(wb, wm)
    dm = comoving_transverse(Ef, h, ok, zs)
    rd = rd_aubourg(wb, wm)
    rs = sound_horizon(Ef, h, wb, zs)   # direct integral: 144.1 Mpc at Planck best fit (CAMB 144.4)
    R = np.sqrt(Om) * (100 * h) * dm / C
    lA = np.pi * dm / rs
    d = np.array([R, lA, wb]) - MEAN
    return float(d @ ICOV @ d), rd, R, lA
