import numpy as np
NAMES, BOUNDS, P0 = ["Om"], [(-0.5, 1.5)], [0.3]   # no dark energy; Om may be negative (Janus is Om = 2 q0)
E = lambda z, p: np.sqrt(np.clip(p[0] * (1 + z) ** 3 + (1 - p[0]) * (1 + z) ** 2, 1e-30, None))
OK = lambda p: 1 - p[0]
