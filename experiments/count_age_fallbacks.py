"""Re-runs AGE-MOEA-II of experiment R1 (same seeds) and counts how often the Euclidean fallback of src/moea.py
replaces pymoo's geodesic distance (review round 5). Checks that every re-run reproduces the stored HV.
Writes results_rev1/age_fallbacks.csv."""
import os
import sys
from multiprocessing import Pool

import pandas as pd

os.environ.setdefault("A1_CONFIG", "config_rev1.json")
os.environ.setdefault("A1_RES", "results_rev1")
sys.path.insert(0, os.path.dirname(__file__))
import run_a1    # noqa: E402  (puts src/ on the path)
import moea      # noqa: E402


def one(args):
    key, r = args
    moea.AGE_FALLBACKS.update(calls=0, fallbacks=0)
    row = run_a1.task((key, "AGE-MOEA-II", r))[0]
    return dict(scenario=key, run=r, HV=row["HV"], **moea.AGE_FALLBACKS)


if __name__ == "__main__":
    keys = list(run_a1.CFG["scenarios_R1"])
    with Pool(4) as pool:
        D = pd.DataFrame(pool.map(one, [(k, r) for k in keys for r in range(run_a1.CFG["runs"])]))
    ref = pd.read_csv(os.path.join(run_a1.RES, "raw_runs.csv"))
    ref = ref[ref.algorithm == "AGE-MOEA-II"][["scenario", "run", "HV"]]
    D = D.merge(ref, on=["scenario", "run"], suffixes=("", "_stored"))
    D["share"] = D.fallbacks / D.calls
    D.to_csv(os.path.join(run_a1.RES, "age_fallbacks.csv"), index=False)
    print("max |HV - stored|:", (D.HV - D.HV_stored).abs().max())
    print("runs with a fallback:", (D.fallbacks > 0).sum(), "of", len(D))
    print("share of distance calls that fell back: overall", D.fallbacks.sum() / D.calls.sum())
    print(D.groupby("scenario")[["calls", "fallbacks", "share"]].mean().round(4).to_string())
