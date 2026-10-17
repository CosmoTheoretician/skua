import numpy as np
NAMES, BOUNDS, P0 = ["Om", "xi"], [(0.05, 0.7), (-1.0, 1.0)], [0.3, 0.0]   # Q = xi H rho_dm into Lambda-like DE, flat
OB = 0.049
def E(z, p):
    om, xi = p
    odm = om - OB; a = 1 / (1 + z); n = 3 - xi
    rdm = odm * a ** (-n)
    rde = (1 - om) - xi * odm * (a ** (-n) - 1) / n
    return np.sqrt(np.clip(OB * (1 + z) ** 3 + rdm + rde, 1e-30, None))
OK = lambda p: 0.0
