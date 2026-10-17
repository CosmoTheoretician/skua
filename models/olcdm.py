import numpy as np
NAMES, BOUNDS, P0 = ["Om", "OL"], [(0.05, 0.7), (0.0, 1.2)], [0.3, 0.7]
E = lambda z, p: np.sqrt(p[0] * (1 + z) ** 3 + (1 - p[0] - p[1]) * (1 + z) ** 2 + p[1])
OK = lambda p: 1 - p[0] - p[1]
