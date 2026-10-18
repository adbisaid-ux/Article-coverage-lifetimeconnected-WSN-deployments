"""Published repair baselines: every repaired random deployment is connected to the sink and stays in the field."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from model import bfs_hops            # noqa: E402
from repairs_pub import dd_mcsa, ky_rh  # noqa: E402

rng = np.random.default_rng(1)
bad = {"KY-RH": 0, "DD-MCSA": 0}
for _ in range(300):
    N, L, rc = rng.choice([30, 50]), rng.choice([100, 140]), rng.choice([8.0, 12.0])
    sink = np.array([L / 2, L / 2])
    xy = rng.uniform(0, L, (N, 2))
    for name, f in [("KY-RH", ky_rh), ("DD-MCSA", dd_mcsa)]:
        y = f(xy, sink, rc, rng)
        hop, _ = bfs_hops(y, sink, rc * rc)
        bad[name] += int((hop < 0).any() or (y < 0).any() or (y > L).any())
print("failures over 300 random deployments:", bad)
assert sum(bad.values()) == 0
