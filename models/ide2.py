import numpy as np
# Interacting dark energy: Q = xi H rho_DM (xi > 0: energy flows DE -> DM), constant w, flat.
OB = 0.049
NAMES, BOUNDS, P0 = ["Om", "w", "xi"], [(0.05, 0.7), (-1.6, -0.5), (-0.3, 0.3)], [0.3, -1.0, 0.0]

def _pieces(p):
    om, w, xi = p
    ode, odm = 1 - om, om - OB
    m = 3 - xi                          # DM dilutes as a^-m
    n = 3 * (1 + w)                     # DE's own dilution
    B = -xi * odm / (n - m) if abs(n - m) > 1e-6 else 0.0
    return ode, odm, m, n, B

def E(z, p):
    ode, odm, m, n, B = _pieces(p)
    a = 1 / (1 + z)
    rdm = odm * a ** -m
    rde = (ode - B) * a ** -n + B * a ** -m
    return np.sqrt(np.clip(OB * a ** -3 + rdm + rde, 1e-30, None))

OK = lambda p: 0.0
def OM(p):                              # the a^-3 matter is baryons only here; DM dilutes as a^-(3-xi).
    ode, odm, m, n, B = _pieces(p)      # For the CMB use the total matter at z ~ 1100 projected to today.
    zs = 1090.0
    return OB + odm * (1 + zs) ** (-xi_of(p)) if False else OB + odm * (1 + zs) ** (m - 3)

def xi_of(p): return p[2]

def VALID(p):
    ode, odm, m, n, B = _pieces(p)
    return (ode - B) > -0.5 and n < 2.9 and abs(p[2]) < 0.3
