"""How far each repair moves sensors (review round 5, item 6). No optimisation runs.
For two scenarios and the first 10 instances of the main campaign, every deployment of the initial population
(100 per instance) is first made connected by the REPAIR control; three sensors chosen at random are then displaced
by N(0, Rc^2) per coordinate, which breaks some links as mutation does, and each repair is applied to the result.
Reported: total displacement / Rc, number of sensors moved, coverage before and after, share fully connected after.
Writes results_pubrate/repair_moves.csv and prints the means."""
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import Problem, agent_step, bfs_hops   # noqa: E402
from repairs_pub import dd_mcsa, ky_rh            # noqa: E402

CFG = json.load(open(os.path.join(ROOT, "experiments", "config.json")))
S = CFG["seeds"]
OUT = os.path.join(ROOT, "results_pubrate")
os.makedirs(OUT, exist_ok=True)
rows = []
for sc in ["N30_Rc8_L140", "N50_Rc12_L100"]:
    s = CFG["scenarios"][sc]
    for r in range(10):
        prob = Problem(s["N"], s["T"], s["L"], CFG["rs"], s["Rc"], target_seed=S["target_seed"][r])
        pop0 = np.random.default_rng(S["init_seed"][r]).uniform(0, s["L"], (CFG["NP"], prob.dim))
        rng = np.random.default_rng(777000 + r)
        for x in pop0:
            rep = lambda xy: agent_step(xy, prob.targets, prob.sink, prob.rs, prob.rc, 0.0, prob.L,   # noqa: E731
                                        prob.rs, rng.random(2 * prob.N), 2)
            base = rep(x.reshape(-1, 2))
            pert = base.copy()
            k = rng.choice(prob.N, 3, replace=False)
            pert[k] = np.clip(pert[k] + rng.normal(0, prob.rc, (3, 2)), 0, prob.L)
            cr0 = prob.evaluate(pert.ravel())[0][0, 0]
            for name, f in [("REPAIR", rep), ("KY-RH", lambda xy: ky_rh(xy, prob.sink, prob.rc, rng)),
                            ("DD-MCSA", lambda xy: dd_mcsa(xy, prob.sink, prob.rc, rng))]:
                y = np.clip(f(pert.copy()), 0, prob.L)
                d = np.sqrt(((y - pert) ** 2).sum(1))
                hop, _ = bfs_hops(y, prob.sink, prob.rc ** 2)
                rows.append(dict(scenario=sc, instance=r, repair=name, disp_over_rc=d.sum() / prob.rc,
                                 moved=int((d > 1e-9).sum()), CR_before=cr0,
                                 CR_after=prob.evaluate(y.ravel())[0][0, 0], all_connected=bool((hop > 0).all())))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(OUT, "repair_moves.csv"), index=False)
print(D.groupby(["scenario", "repair"])[["disp_over_rc", "moved", "CR_before", "CR_after", "all_connected"]]
      .mean().round(3).to_string())
