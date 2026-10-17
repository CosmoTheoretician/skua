import numpy as np
NAMES, BOUNDS, P0 = ["Om"], [(0.05, 0.7)], [0.3]
E = lambda z, p: np.sqrt(p[0] * (1 + z) ** 3 + 1 - p[0])
OK = lambda p: 0.0
