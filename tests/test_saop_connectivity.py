"""Property check for SA-op: after one application every sensor is connected to the sink."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from model import Problem, agent_step, bfs_hops, exact_c, evaluate_batch   # noqa: E402


def test_saop_connectivity(n=300):
    for s in range(n):
        rng = np.random.default_rng(s)
        N, L, rc = [30, 50][s % 2], [100, 140][(s // 2) % 2], [8, 12][(s // 4) % 2]
        P = Problem(N, 100, L, 8, rc, target_seed=s)
        y = agent_step(rng.uniform(0, L, (N, 2)), P.targets, P.sink, 8.0, float(rc), 0.0, float(L), 8.0,
                       rng.random(2 * N))
        h, _ = bfs_hops(y, P.sink, rc * rc)
        assert (h > 0).all(), s


def test_ctree_ge_cstar(n=50):
    for s in range(n):
        rng = np.random.default_rng(1000 + s)
        P = Problem(50, 100, 100, 8, 12, target_seed=s)
        x = rng.uniform(30, 70, (1, 100))
        _, c_tree, m = evaluate_batch(x, P.targets, P.sink, 8.0, 12.0)
        if m[0] > 0:
            assert c_tree[0] >= exact_c(x[0].reshape(-1, 2), P.sink, 12.0) - 1e-3


if __name__ == "__main__":
    test_saop_connectivity()
    test_ctree_ge_cstar()
    print("ok")
