import numpy as np
NAMES, BOUNDS, P0 = ["Om"], [(0.05, 0.9)], [0.3]   # flat self-accelerating DGP: Orc = ((1-Om)/2)^2
def E(z, p):
    om = p[0]; orc = ((1 - om) / 2) ** 2
    return np.sqrt(om * (1 + z) ** 3 + orc) + np.sqrt(orc)
OK = lambda p: 0.0
