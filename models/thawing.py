import numpy as np
NAMES, BOUNDS, P0 = ["Om", "w0", "wa"], [(0.05, 0.7), (-1.0, 0.0), (-3.0, 3.0)], [0.3, -0.9, 0.0]
def E(z, p):
    om, w0, wa = p
    if w0 + wa < -1:
        return np.full_like(np.atleast_1d(z), np.nan, dtype=float)
    f = (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))
    return np.sqrt(om * (1 + z) ** 3 + (1 - om) * f)
OK = lambda p: 0.0
