import numpy as np
NAMES = ["Om", "f05", "f1", "f2"]                         # dark-energy density relative to today at z = 0.5, 1, 2
BOUNDS, P0 = [(0.05, 0.7), (-3, 4), (-3, 4), (-3, 4)], [0.3, 1, 1, 1]
ZN = np.array([0.0, 0.5, 1.0, 2.0, 2.4])
def E(z, p):
    om = p[0]; f = np.interp(z, ZN, [1.0, p[1], p[2], p[3], p[3]])
    return np.sqrt(np.clip(om * (1 + z) ** 3 + (1 - om) * f, 1e-30, None))
OK = lambda p: 0.0
