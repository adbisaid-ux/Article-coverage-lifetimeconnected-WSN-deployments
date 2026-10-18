"""Revision runs: statistics for R1, R2 (SA-op sensitivity), R3 (EPS-DE+REPAIR), L (long runs, equal time)
and the evaluation timing quoted in the paper. Paired Wilcoxon signed-rank tests, Holm correction within
each family stated below.
Writes results_rev1/stats_R1.csv, sensitivity_saop.csv, stats_R3.csv, results_long/stats_L.csv,
equal_time.csv, growth_L.csv and results_rev1/timing_eval.json."""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import Problem   # noqa: E402

M = pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv"))
S = pd.read_csv(os.path.join(ROOT, "results_rev1", "raw_runs.csv"))
L = pd.read_csv(os.path.join(ROOT, "results_long", "raw_runs.csv"))
H = pd.read_csv(os.path.join(ROOT, "results_long", "hv_history.csv"))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(run, 1.0)
    return adj


def verdict(r):
    return np.where(r.p_holm < 0.05, np.where(r.ma > r.mb, "+", "-"), "=")


# R1: NSGA-II+SA vs AGE-MOEA-II and EPS-DE; pymoo NSGA-II with vs without duplicate elimination.
# Holm family: the 24 tests of one metric.
x = pd.concat([M, S])
rows = []
for k, d in x.groupby("scenario"):
    for met in ["HV", "max_CR"]:
        w = d.pivot(index="run", columns="algorithm", values=met)
        for o in ["AGE-MOEA-II", "EPS-DE", "NSGA-II(pymoo,nodup)"]:
            a, b = ("NSGA-II+SA", o) if o != "NSGA-II(pymoo,nodup)" else ("NSGA-II(pymoo,nodup)", "NSGA-II(pymoo)")
            rows.append(dict(scenario=k, metric=met, a=a, b=b, ma=w[a].mean(), mb=w[b].mean(),
                             p=stats.wilcoxon(w[a], w[b]).pvalue))
r = pd.DataFrame(rows)
for met in ["HV", "max_CR"]:
    q = r.metric == met
    r.loc[q, "p_holm"] = holm(r.loc[q, "p"])
r["res"] = verdict(r)
r["gain%"] = 100 * (r.ma / r.mb - 1)
r.to_csv(os.path.join(ROOT, "results_rev1", "stats_R1.csv"), index=False)

# R2: sensitivity of NSGA-II+SA to the step delta and the projection margin. Holm family: 24 tests.
x = pd.concat([M[M.algorithm == "NSGA-II+SA"], S[S.algorithm.str.startswith("NSGA-II+SA ")]])
rows = []
for k, d in x.groupby("scenario"):
    w = d.pivot(index="run", columns="algorithm", values="HV")
    for o in ["NSGA-II+SA d=4", "NSGA-II+SA d=16", "NSGA-II+SA m=0.8"]:
        rows.append(dict(scenario=k, other=o, base=w["NSGA-II+SA"].mean(), var=w[o].mean(),
                         p=stats.wilcoxon(w["NSGA-II+SA"], w[o]).pvalue))
r = pd.DataFrame(rows)
r["p_holm"] = holm(r.p)
r["diff_pct"] = 100 * (r["var"] / r.base - 1)
r.to_csv(os.path.join(ROOT, "results_rev1", "sensitivity_saop.csv"), index=False)

# R3: EPS-DE+REPAIR vs EPS-DE, NSGA-II+SA and NSGA-II+REPAIR (main and E1 campaigns). Holm family: 24 tests per metric.
B = pd.read_csv(os.path.join(ROOT, "results_boucle", "raw_runs.csv"))
x = pd.concat([M, S, B[B.algorithm == "NSGA-II+REPAIR"]])
if (x.algorithm == "EPS-DE+REPAIR").sum() == 240:
    rows = []
    for k, d in x.groupby("scenario"):
        if "EPS-DE+REPAIR" not in set(d.algorithm):
            continue
        for met in ["HV", "max_CR"]:
            w = d.pivot(index="run", columns="algorithm", values=met)
            for a, b in [("EPS-DE+REPAIR", "EPS-DE"), ("NSGA-II+SA", "EPS-DE+REPAIR"),
                         ("NSGA-II+REPAIR", "EPS-DE+REPAIR")]:
                rows.append(dict(scenario=k, metric=met, a=a, b=b, ma=w[a].mean(), mb=w[b].mean(),
                                 p=stats.wilcoxon(w[a], w[b]).pvalue))
    r = pd.DataFrame(rows)
    for met in ["HV", "max_CR"]:
        q = r.metric == met
        r.loc[q, "p_holm"] = holm(r.loc[q, "p"])
    r["res"] = verdict(r)
    r["gain%"] = 100 * (r.ma / r.mb - 1)
    r.to_csv(os.path.join(ROOT, "results_rev1", "stats_R3.csv"), index=False)

# L: 5e5 evaluations. Holm family: 16 tests.
rows = []
for k, d in L.groupby("scenario"):
    w = d.pivot(index="run", columns="algorithm", values="HV")
    for o in ["NSGA-II", "NSGA-II+REPAIR"]:
        rows.append(dict(scenario=k, b=o, ma=w["NSGA-II+SA"].mean(), mb=w[o].mean(),
                         p=stats.wilcoxon(w["NSGA-II+SA"], w[o]).pvalue))
r = pd.DataFrame(rows)
r["p_holm"] = holm(r.p)
r["gain%"] = 100 * (r.ma / r.mb - 1)
r["res"] = verdict(r)
r.to_csv(os.path.join(ROOT, "results_long", "stats_L.csv"), index=False)


def hv_at(fe):
    return H[H.fe <= fe].sort_values("fe").groupby(["scenario", "algorithm", "run"]).HV.last().reset_index()


# Equal time: NSGA-II+SA at 1e5 evaluations vs NSGA-II at 1.3e5 (mean run time ratio about 1.3).
a1, a13 = hv_at(100000), hv_at(130000)
rows = []
for k in sorted(L.scenario.unique()):
    a = a1[(a1.scenario == k) & (a1.algorithm == "NSGA-II+SA")].set_index("run").HV
    b = a13[(a13.scenario == k) & (a13.algorithm == "NSGA-II")].set_index("run").HV
    rows.append(dict(scenario=k, SA_100k=a.mean(), NSGA2_130k=b.mean(), p=stats.wilcoxon(a, b).pvalue))
pd.DataFrame(rows).to_csv(os.path.join(ROOT, "results_long", "equal_time.csv"), index=False)

# Growth of mean HV from 1e5 to 5e5 evaluations: per scenario (growth_L.csv) and for the mean over all
# 240 runs of each algorithm (the figure quoted in the paper).
g = pd.merge(a1.groupby(["scenario", "algorithm"]).HV.mean(), hv_at(500000).groupby(["scenario", "algorithm"]).HV.mean(),
             left_index=True, right_index=True, suffixes=("_1e5", "_5e5")).reset_index()
g["growth%"] = 100 * (g.HV_5e5 / g.HV_1e5 - 1)
g.to_csv(os.path.join(ROOT, "results_long", "growth_L.csv"), index=False)
pooled = g.groupby("algorithm")[["HV_1e5", "HV_5e5"]].sum()
print("growth of mean HV over all runs (%):", (100 * (pooled.HV_5e5 / pooled.HV_1e5 - 1)).round(1).to_dict())

# Evaluation timing: one call of Problem.evaluate on a batch of 100 random deployments (N50, Rc12, L100).
prob = Problem(50, 100, 100, 10, 12, target_seed=20261016)
pop = np.random.default_rng(0).uniform(0, 100, (100, prob.dim))
prob.evaluate(pop)   # JIT warm-up
t0 = time.perf_counter()
for _ in range(200):
    prob.evaluate(pop)
us = 1e6 * (time.perf_counter() - t0) / (200 * 100)
json.dump(dict(scenario="N50_Rc12_L100", batch=100, repeats=200, us_per_eval=us),
          open(os.path.join(ROOT, "results_rev1", "timing_eval.json"), "w"), indent=1)
print(f"evaluation: {us:.1f} us per deployment")
