import numpy as np
NAMES, BOUNDS, P0 = ["q0"], [(-0.21, 0.6)], [0.0]
E = lambda z, p: (1 + z) * np.sqrt(np.clip(1 + 2 * p[0] * z, 1e-15, None))
OK = lambda p: 1 - 2 * p[0]
