import numpy as np
NAMES, BOUNDS, P0 = ["Om", "Ox", "w"], [(0.0, 1.0), (-1.5, 1.5), (-2.5, 0.5)], [0.3, 0.7, -1.0]
def E(z, p):
    om, ox, w = p
    return np.sqrt(np.clip(om * (1 + z) ** 3 + (1 - om - ox) * (1 + z) ** 2 + ox * (1 + z) ** (3 * (1 + w)), 1e-30, None))
OK = lambda p: 1 - p[0] - p[1]
