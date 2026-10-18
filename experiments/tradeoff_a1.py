"""C3: operational reading of the NSGA-II+SA fronts. For each scenario and each relay-load
level c (lifetime T(c)), the largest connected coverage reached on the final front of each
run; median and IQR over 30 runs. Writes results/tables/T7_tradeoff.{csv,md}."""
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import lifetime_from_c, radio_costs   # noqa: E402

CFG = json.load(open(os.path.join(ROOT, "experiments", "config.json")))
fr = pd.read_csv(os.path.join(ROOT, "results", "fronts.csv"))
fr = fr[fr.algorithm == "NSGA-II+SA"]
rows = []
for k, s in CFG["scenarios"].items():
    e_tx, e_rx = radio_costs(s["Rc"])
    Tmax = lifetime_from_c(1.0, s["Rc"])
    d = fr[fr.scenario == k].copy()
    T = d.Tn * Tmax
    d["c"] = np.round((0.5 / T + e_rx) / (e_tx + e_rx), 3)
    for c in [1, 2, 3, 4, 5, 6, 8, 10]:
        best = d[d.c <= c + 1e-6].groupby("run").CR.max().reindex(range(CFG["runs"])).fillna(0.0)
        rows.append(dict(scenario=k, c=c, T_rounds=int(lifetime_from_c(c, s["Rc"])),
                         CR_median=100 * best.median(), CR_q1=100 * best.quantile(.25),
                         CR_q3=100 * best.quantile(.75)))
t = pd.DataFrame(rows)
t.to_csv(os.path.join(ROOT, "results", "tables", "T7_tradeoff.csv"), index=False)
piv = t.assign(v=t.apply(lambda r: f"{r.CR_median:.0f} [{r.CR_q1:.0f}, {r.CR_q3:.0f}]", axis=1)) \
       .pivot(index="scenario", columns="c", values="v")[[1, 2, 3, 4, 5, 6, 8, 10]]
rounds = t[t.scenario == "N50_Rc12_L100"].set_index("c").T_rounds
piv.columns = [f"c={c} ({rounds[c]} r)" for c in piv.columns]
open(os.path.join(ROOT, "results", "tables", "T7_tradeoff.md"), "w").write(
    "### Table 7. NSGA-II+SA: largest connected coverage (%) reachable with maximum relay load c "
    "(lifetime in rounds for Rc = 12 m), median [Q1, Q3] over 30 runs\n\n" + piv.to_markdown() +
    "\n\nRounds for Rc = 8 m: " + ", ".join(f"c={c}: {int(lifetime_from_c(c, 8))}" for c in [1, 2, 3, 4, 5, 6, 8, 10]) + ".\n")
print(open(os.path.join(ROOT, "results", "tables", "T7_tradeoff.md")).read())
# distinct c values present
print(sorted(set(np.round((0.5 / (fr.Tn * 2430.0155520995336) + 2e-4) / (2.0576e-4 + 2e-4), 2)))[:15])
