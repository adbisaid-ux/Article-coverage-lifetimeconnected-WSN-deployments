"""Analysis of the extra campaigns of review round 4.
V  (results_var): intra-instance variability. 3 instances x 10 algorithm seeds per scenario.
TU (results_tune): NSGA-II variation-parameter grid on 2 scenarios, 10 runs with disjoint seeds.
TN (results_tuned, if present): the selected NSGA-II setting on the main campaign seeds, against NSGA-II+SA.
Writes the tables into each results folder and prints a summary."""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(run, 1.0)
    return adj


# ---------------------------------------------------------------- V: variability
V = pd.read_csv(os.path.join(ROOT, "results_var", "raw_runs.csv"))
V["instance"] = V.run // 10
M = pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv"))
B = pd.read_csv(os.path.join(ROOT, "results_boucle", "raw_runs.csv"))
rows = []
for (sc, alg), d in V.groupby(["scenario", "algorithm"]):
    g = d.groupby("instance").HV
    within = np.sqrt(g.var(ddof=1).mean())                 # pooled within-instance SD (algorithm seeds)
    between = g.mean().std(ddof=1)                         # SD of the 3 instance means
    ref = M if alg != "NSGA-II+REPAIR" else B
    across = ref[(ref.scenario == sc) & (ref.algorithm == alg)].HV.std(ddof=1)   # 30 instances, main seeds
    rows.append(dict(scenario=sc, algorithm=alg, sd_within_instance=within, sd_between_3_instances=between,
                     sd_main_30_instances=across, share_within=within ** 2 / across ** 2 if across > 0 else np.nan))
VV = pd.DataFrame(rows)
VV.to_csv(os.path.join(ROOT, "results_var", "variance.csv"), index=False)

rows = []
for (sc, inst), d in V.groupby(["scenario", "instance"]):
    w = d.pivot(index="run", columns="algorithm", values="HV")
    for o in ["NSGA-II", "NSGA-II+REPAIR"]:
        diff = w["NSGA-II+SA"] - w[o]
        p = stats.wilcoxon(w["NSGA-II+SA"], w[o]).pvalue if (diff != 0).any() else 1.0
        rows.append(dict(scenario=sc, instance=inst, other=o, mean_diff=diff.mean(),
                         gain_pct=100 * (w["NSGA-II+SA"].mean() / w[o].mean() - 1),
                         seeds_SA_better=int((diff > 0).sum()), p=p))
VI = pd.DataFrame(rows)
for o in ["NSGA-II", "NSGA-II+REPAIR"]:
    q = VI.other == o
    VI.loc[q, "p_holm"] = holm(VI.loc[q, "p"])             # Holm over the 24 instances of each comparison
VI["res"] = np.where(VI.p_holm < 0.05, np.where(VI.mean_diff > 0, "+", "-"), "=")
VI.to_csv(os.path.join(ROOT, "results_var", "per_instance.csv"), index=False)
print("V: per-instance verdicts (Holm over 24 instances)")
print(VI.groupby(["other", "res"]).size().unstack(fill_value=0))
print("V: SA better on", VI[VI.other == "NSGA-II"].seeds_SA_better.sum(), "of 240 seed pairs vs NSGA-II;",
      VI[VI.other == "NSGA-II+REPAIR"].seeds_SA_better.sum(), "of 240 vs REPAIR")
print(VV.round(4).to_string())

# ---------------------------------------------------------------- TU: tuning grid
p = os.path.join(ROOT, "results_tune", "raw_runs.csv")
if os.path.exists(p):
    T = pd.read_csv(p)
    G = T.groupby("algorithm").HV.agg(["mean", "std", "count"]).sort_values("mean", ascending=False)
    G["mean_by_scenario"] = [T[T.algorithm == a].groupby("scenario").HV.mean().round(4).to_dict() for a in G.index]
    G.to_csv(os.path.join(ROOT, "results_tune", "grid.csv"))
    best = G.index[0]
    json.dump(dict(rule="highest mean HV over the 20 tuning runs", selected=best,
                   default="NSGA-II ec=15 em=20", n_runs=int(G.loc[best, "count"])),
              open(os.path.join(ROOT, "results_tune", "selected.json"), "w"), indent=1)
    print("TU grid:\n", G.round(4).to_string(), "\nselected:", best)

# ---------------------------------------------------------------- TN: tuned NSGA-II vs NSGA-II+SA (main seeds)
p = os.path.join(ROOT, "results_tuned", "raw_runs.csv")
if os.path.exists(p):
    TN = pd.read_csv(p)
    alg = TN.algorithm.iloc[0]
    x = pd.concat([M[M.algorithm.isin(["NSGA-II+SA", "NSGA-II"])], TN])
    rows = []
    for sc, d in x.groupby("scenario"):
        for met in ["HV", "max_CR"]:
            w = d.pivot(index="run", columns="algorithm", values=met)
            for a, b in [("NSGA-II+SA", alg), (alg, "NSGA-II")]:
                rows.append(dict(scenario=sc, metric=met, a=a, b=b, ma=w[a].mean(), mb=w[b].mean(),
                                 p=stats.wilcoxon(w[a], w[b]).pvalue))
    R = pd.DataFrame(rows)
    for (met, a), q in R.groupby(["metric", "a"]).groups.items():
        R.loc[q, "p_holm"] = holm(R.loc[q, "p"])
    R["res"] = np.where(R.p_holm < 0.05, np.where(R.ma > R.mb, "+", "-"), "=")
    R["gain%"] = 100 * (R.ma / R.mb - 1)
    R.to_csv(os.path.join(ROOT, "results_tuned", "stats_tuned.csv"), index=False)
    print("TN:\n", R.round(4).to_string())

# ---------------------------------------------------------------- PB / PA: published repairs
for folder, base, regime in [("results_pub", None, "raw"), ("results_pub_agg", "results_agg", "agg")]:
    p = os.path.join(ROOT, folder, "raw_runs.csv")
    if not os.path.exists(p):
        continue
    P = pd.read_csv(p)
    if regime == "raw":
        ref = pd.concat([M[M.algorithm.isin(["NSGA-II+SA", "NSGA-II"])], B[B.algorithm == "NSGA-II+REPAIR"]])
    else:
        A = pd.read_csv(os.path.join(ROOT, base, "raw_runs.csv"))
        ref = A[A.scenario.str.startswith("A_") & A.algorithm.isin(["NSGA-II+SA", "NSGA-II", "NSGA-II+REPAIR"])]
    x = pd.concat([ref, P])
    rows = []
    for sc, d in x.groupby("scenario"):
        for met in ["HV", "max_CR"]:
            w = d.pivot(index="run", columns="algorithm", values=met)
            for a, b in [("NSGA-II+SA", "NSGA-II+KY-RH"), ("NSGA-II+SA", "NSGA-II+DD-MCSA"),
                         ("NSGA-II+REPAIR", "NSGA-II+KY-RH"), ("NSGA-II+REPAIR", "NSGA-II+DD-MCSA"),
                         ("NSGA-II+KY-RH", "NSGA-II"), ("NSGA-II+DD-MCSA", "NSGA-II")]:
                pv = stats.wilcoxon(w[a], w[b]).pvalue if (w[a] != w[b]).any() else 1.0
                rows.append(dict(scenario=sc, metric=met, a=a, b=b, ma=w[a].mean(), mb=w[b].mean(), p=pv,
                                 a12=(np.sign(w[a].values[:, None] - w[b].values[None, :]) + 1).mean() / 2))
    R = pd.DataFrame(rows)
    for met in ["HV", "max_CR"]:
        q = R.metric == met
        R.loc[q, "p_holm"] = holm(R.loc[q, "p"])          # Holm over the 48 tests of one metric
    R["res"] = np.where(R.p_holm < 0.05, np.where(R.ma > R.mb, "+", "-"), "=")
    R["gain%"] = 100 * (R.ma / R.mb - 1)
    R.to_csv(os.path.join(ROOT, folder, "stats_pub.csv"), index=False)
    print(folder)
    print(R.groupby(["metric", "a", "b", "res"]).size().unstack(fill_value=0))
    print(R[R.metric == "HV"].groupby(["a", "b"])["gain%"].agg(["min", "max"]).round(1))
    print(P.groupby("algorithm").time_s.mean().round(1))

# ---------------------------------------------------------------- PR: repair rate sensitivity (review round 5)
p = os.path.join(ROOT, "results_pubrate", "raw_runs.csv")
if os.path.exists(p):
    P = pd.read_csv(p)
    Q = pd.read_csv(os.path.join(ROOT, "results_pub", "raw_runs.csv"))
    x = pd.concat([M[M.algorithm.isin(["NSGA-II+SA", "NSGA-II+SA p_sa=0.5", "NSGA-II"])],
                   B[B.algorithm == "NSGA-II+REPAIR"], Q, P])
    x = x[x.scenario.isin(P.scenario.unique())]
    x.groupby(["scenario", "algorithm"])[["HV", "max_CR", "time_s"]].mean().to_csv(
        os.path.join(ROOT, "results_pubrate", "means.csv"))
    rows = []
    for sc, d in x.groupby("scenario"):
        for met in ["HV", "max_CR"]:
            w = d.pivot(index="run", columns="algorithm", values=met)
            for b in [c for c in w.columns if "KY-RH" in c or "DD-MCSA" in c]:
                rows.append(dict(scenario=sc, metric=met, a="NSGA-II+SA", b=b, ma=w["NSGA-II+SA"].mean(),
                                 mb=w[b].mean(), p=stats.wilcoxon(w["NSGA-II+SA"], w[b]).pvalue))
    R = pd.DataFrame(rows)
    for met in ["HV", "max_CR"]:
        q = R.metric == met
        R.loc[q, "p_holm"] = holm(R.loc[q, "p"])          # Holm over the 12 tests of one metric
    R["res"] = np.where(R.p_holm < 0.05, np.where(R.ma > R.mb, "+", "-"), "=")
    R["gain%"] = 100 * (R.ma / R.mb - 1)
    R.to_csv(os.path.join(ROOT, "results_pubrate", "stats_rate.csv"), index=False)
    print("PR:\n", R.round(4).to_string())
