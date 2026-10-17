import numpy as np
NAMES, BOUNDS, P0 = ["Om", "w"], [(0.05, 0.7), (-2.5, 0.0)], [0.3, -1.0]
E = lambda z, p: np.sqrt(p[0] * (1 + z) ** 3 + (1 - p[0]) * (1 + z) ** (3 * (1 + p[1])))
OK = lambda p: 0.0
