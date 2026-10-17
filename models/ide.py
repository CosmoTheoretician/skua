import numpy as np
# Interacting dark energy: Q = xi H rho_DE (xi > 0: energy flows DE -> DM), constant w, flat.
OB = 0.049
NAMES, BOUNDS, P0 = ["Om", "w", "xi"], [(0.05, 0.7), (-1.6, -0.5), (-0.3, 0.3)], [0.3, -1.0, 0.0]

def _pieces(p):
    om, w, xi = p
    ode, odm = 1 - om, om - OB
    n = 3 * (1 + w) + xi
    A = xi * ode / (3 - n) if abs(3 - n) > 1e-6 else 0.0
    return ode, odm, n, A

def E(z, p):
    ode, odm, n, A = _pieces(p)
    a = 1 / (1 + z)
    rdm = (odm - A) * a ** -3 + A * a ** -n
    return np.sqrt(np.clip(OB * a ** -3 + rdm + ode * a ** -n, 1e-30, None))

OK = lambda p: 0.0
def OM(p):                     # matter density that scales as a^-3 (what the CMB and r_d see)
    ode, odm, n, A = _pieces(p)
    return OB + odm - A

def VALID(p):                  # early dark matter must be positive and dilute faster than the DE term
    ode, odm, n, A = _pieces(p)
    return (odm - A) > 0 and n < 2.9
