"""Revision round 1, analyses on existing data (main campaign + /boucle):
 - distinct objective vectors per final front;
 - IGD+ w.r.t. the per-instance reference front (non-dominated union of all algorithms' fronts);
 - HV without the no-relay point (c = 1), i.e. the part of the front that uses relays;
 - largest coverage per load level c for all algorithms (median over runs);
 - paired differences NSGA-II+SA minus each competitor with 95 % CI, and share of HV variance
   explained by the instance (run) vs the algorithm;
 - max_CR tests among control operators (SA, RND, REPAIR) with Holm over the 24 tests;
 - fig1 redrawn on a single run (same instance for all algorithms).
Writes results/tables_rev1/*.md|csv and figures/fig1_fronts_same_run.{pdf,png}."""
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import lifetime_from_c, radio_costs   # noqa: E402
from moea import hv2d, nondominated               # noqa: E402

CFG = json.load(open(os.path.join(ROOT, "experiments", "config.json")))
OUT = os.path.join(ROOT, "results", "tables_rev1")
os.makedirs(OUT, exist_ok=True)
K8 = list(CFG["scenarios"])
MAIN = CFG["algorithms"]
fr = pd.concat([pd.read_csv(os.path.join(ROOT, "results", "fronts.csv")),
                pd.read_csv(os.path.join(ROOT, "results_boucle", "fronts.csv"))])
raw = pd.concat([pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv")),
                 pd.read_csv(os.path.join(ROOT, "results_boucle", "raw_runs.csv"))])
fr = fr[fr.scenario.isin(K8)]
ALL = [a for a in fr.algorithm.unique()]


def load_c(scen, Tn):
    rc = CFG["scenarios"][scen]["Rc"]
    e_tx, e_rx = radio_costs(rc)
    T = Tn * lifetime_from_c(1.0, rc)
    return np.round((CFG["radio"]["E0_J"] / T + e_rx) / (e_tx + e_rx)).astype(int)


def igd_plus(A, R):
    """IGD+ (minimisation) of approximation A w.r.t. reference set R."""
    d = np.maximum(A[None, :, :] - R[:, None, :], 0.0)
    return np.sqrt((d ** 2).sum(-1)).min(1).mean()


rows = []
for (k, r), d in fr.groupby(["scenario", "run"]):
    G_all = 1.0 - d[["CR", "Tn"]].values
    ref = np.unique(G_all[nondominated(G_all)], axis=0)
    for a, da in d.groupby("algorithm"):
        G = 1.0 - da[["CR", "Tn"]].values
        Gd = np.unique(G, axis=0)
        c = load_c(k, da.Tn.values)
        relay = G[c >= 2]
        rows.append(dict(scenario=k, run=r, algorithm=a, n_distinct=len(Gd), IGDp=igd_plus(Gd, ref),
                         HV=hv2d(G), HV_relay=hv2d(relay) if len(relay) else 0.0,
                         **{f"CR_c{cc}": (da.CR[c <= cc].max() if (c <= cc).any() else 0.0)
                            for cc in [1, 2, 4, 6, 10]}))
M = pd.DataFrame(rows)
M["share_c1"] = 1 - M.HV_relay / M.HV
M.to_csv(os.path.join(OUT, "front_metrics.csv"), index=False)

order = MAIN + [a for a in ALL if a not in MAIN]
summ = M.groupby("algorithm").agg(n_distinct=("n_distinct", "mean"), IGDp=("IGDp", "mean"),
                                  HV_relay=("HV_relay", "mean"), share_c1=("share_c1", "mean")).reindex(order)
open(os.path.join(OUT, "T8_front_metrics_overall.md"), "w").write(
    "### Front metrics averaged over the 8 scenarios x 30 runs\n\n" + summ.round(4).to_markdown() +
    "\n\nn_distinct: distinct objective vectors on the final front; IGD+: w.r.t. the non-dominated union of all "
    "algorithms on the same instance (lower is better); HV_relay: HV of the front points with load c >= 2; "
    "share_c1: share of HV contributed by the no-relay point.\n")


def per_scen_table(metric, fmt, algs):
    t = M[M.algorithm.isin(algs)].groupby(["scenario", "algorithm"])[metric].agg(["mean", "std"])
    tab = pd.DataFrame({a: [f"{t.loc[(k, a), 'mean']:{fmt}} ± {t.loc[(k, a), 'std']:{fmt}}" for k in K8]
                        for a in algs}, index=K8)
    tab.index.name = "Scenario"
    return tab


# statistical tests on IGD+ (lower better) and HV_relay for the main algorithms
def tests(metric, better_low):
    out = []
    for k in K8:
        w = M[(M.scenario == k) & M.algorithm.isin(MAIN)].pivot(index="run", columns="algorithm", values=metric)
        for o in MAIN[1:]:
            x, y = w["NSGA-II+SA"], w[o]
            p = 1.0 if (x == y).all() else stats.wilcoxon(x, y).pvalue
            out.append(dict(scenario=k, other=o, p=p, mx=x.mean(), my=y.mean()))
    st = pd.DataFrame(out)
    o = np.argsort(st.p.values)
    adj = np.empty(len(st))
    run = 0
    for i, j in enumerate(o):
        run = max(run, (len(st) - i) * st.p.values[j])
        adj[j] = min(run, 1)
    st["p_holm"] = adj
    good = (st.mx < st.my) if better_low else (st.mx > st.my)
    st["result"] = np.where(st.p_holm < 0.05, np.where(good, "+", "-"), "=")
    return st


for metric, fmt, low in [("IGDp", ".4f", True), ("HV_relay", ".4f", False)]:
    st = tests(metric, low)
    st.to_csv(os.path.join(OUT, f"stats_{metric}.csv"), index=False)
    tab = per_scen_table(metric, fmt, MAIN)
    for k in K8:
        for o in MAIN[1:]:
            tab.loc[k, o] += f" ({st[(st.scenario == k) & (st.other == o)].result.iloc[0]})"
    open(os.path.join(OUT, f"T9_{metric}.md"), "w").write(
        f"### {metric} per scenario, mean ± sd over 30 runs\n\n" + tab.to_markdown() +
        f"\n\n(+/−/=): NSGA-II+SA significantly better / worse / tie (paired Wilcoxon, Holm over {len(st)}).\n"
        f"Summary: {''.join(st.result)}\n")

# coverage per load level, all main algorithms (median)
cl = M[M.algorithm.isin(MAIN)].groupby(["scenario", "algorithm"])[["CR_c2", "CR_c4", "CR_c6", "CR_c10"]].median()
cl = (100 * cl).round(0).astype(int).astype(str)
cl = cl.CR_c2 + " / " + cl.CR_c4 + " / " + cl.CR_c6 + " / " + cl.CR_c10
tab = cl.unstack()[MAIN]
tab.index.name = "Scenario"
open(os.path.join(OUT, "T10_coverage_per_level.md"), "w").write(
    "### Median largest coverage (%) reachable with load c = 2 / 4 / 6 / 10, all algorithms\n\n" + tab.to_markdown() + "\n")

# paired differences with 95% CI and variance decomposition
rows = []
for k in K8:
    w = raw[(raw.scenario == k) & raw.algorithm.isin(MAIN)].pivot(index="run", columns="algorithm", values="HV")[MAIN]
    for o in MAIN[1:]:
        dd = w["NSGA-II+SA"] - w[o]
        h = stats.t.ppf(0.975, len(dd) - 1) * dd.std(ddof=1) / np.sqrt(len(dd))
        rows.append(dict(scenario=k, other=o, mean_diff=dd.mean(), ci_low=dd.mean() - h, ci_high=dd.mean() + h))
    # two-way decomposition (run = instance, algorithm), no interaction replicate
    grand = w.values.mean()
    ss_tot = ((w.values - grand) ** 2).sum()
    ss_run = len(MAIN) * ((w.mean(1) - grand) ** 2).sum()
    ss_alg = len(w) * ((w.mean(0) - grand) ** 2).sum()
    rows.append(dict(scenario=k, other="_variance_", mean_diff=ss_run / ss_tot, ci_low=ss_alg / ss_tot,
                     ci_high=1 - (ss_run + ss_alg) / ss_tot))
P = pd.DataFrame(rows)
P.to_csv(os.path.join(OUT, "paired_diff_ci.csv"), index=False)
V = P[P.other == "_variance_"].set_index("scenario")[["mean_diff", "ci_low", "ci_high"]]
V.columns = ["share_instance", "share_algorithm", "share_residual"]
D = P[(P.other != "_variance_")]
D2 = D.assign(v=D.apply(lambda r: f"{r.mean_diff:.4f} [{r.ci_low:.4f}, {r.ci_high:.4f}]", axis=1)) \
      .pivot(index="scenario", columns="other", values="v")[MAIN[1:]].reindex(K8)
open(os.path.join(OUT, "T11_paired_differences.md"), "w").write(
    "### HV difference NSGA-II+SA minus competitor, mean [95% CI] over 30 paired runs\n\n" + D2.to_markdown() +
    "\n\n### Share of HV sum of squares explained by instance and algorithm (main 7 algorithms)\n\n" +
    V.round(3).to_markdown() + "\n")

# control operators: max_CR tests with Holm over all 24
pairs = [("NSGA-II+SA", "NSGA-II+RND"), ("NSGA-II+SA", "NSGA-II+REPAIR"), ("NSGA-II+RND", "NSGA-II+REPAIR")]
out = []
for k in K8:
    w = raw[raw.scenario == k].pivot(index="run", columns="algorithm", values="max_CR")
    for a, b in pairs:
        out.append(dict(scenario=k, a=a, b=b, ma=w[a].mean(), mb=w[b].mean(), p=stats.wilcoxon(w[a], w[b]).pvalue))
C = pd.DataFrame(out)
o = np.argsort(C.p.values)
adj = np.empty(len(C))
run = 0
for i, j in enumerate(o):
    run = max(run, (len(C) - i) * C.p.values[j])
    adj[j] = min(run, 1)
C["p_holm"] = adj
C["result"] = np.where(C.p_holm < 0.05, np.where(C.ma > C.mb, "+", "-"), "=")
C.to_csv(os.path.join(OUT, "stats_controls_maxCR.csv"), index=False)
lines = [f"{a} vs {b}: {''.join(C[(C.a == a) & (C.b == b)].result)}" for a, b in pairs]
open(os.path.join(OUT, "T12_controls_maxCR.md"), "w").write(
    "### Largest coverage among control operators (paired Wilcoxon, Holm over 24 tests), scenarios in config order\n\n"
    + C.round(4).to_markdown(index=False) + "\n\n" + "\n".join(lines) + "\n")

# fig1 on a single run per scenario: the run with median HV of NSGA-II+SA
cols = {"NSGA-II+SA": "#1f77b4", "NSGA-II": "#ff7f0e", "NSGA-II(pymoo)": "#2ca02c", "NSGA-III": "#d62728",
        "MOEA/D": "#9467bd", "SMS-EMOA": "#8c564b", "DE-CR": "#7f7f7f"}
fig, axes = plt.subplots(2, 4, figsize=(18, 8), sharey=True)
for ax, k in zip(axes.ravel(), K8):
    h = raw[(raw.scenario == k) & (raw.algorithm == "NSGA-II+SA")].set_index("run").HV
    r = (h - h.median()).abs().idxmin()
    for a in MAIN:
        d = fr[(fr.scenario == k) & (fr.algorithm == a) & (fr.run == r)].drop_duplicates(["CR", "Tn"]).sort_values("CR")
        ax.step(100 * d.CR, d.Tn, where="pre", marker="o", ms=3, label=a, color=cols[a], lw=1.2)
    ax.set_title(f"{k} (run {r})", fontsize=10)
    ax.set_xlabel("Connected coverage (%)")
    ax.grid(alpha=0.3)
for ax in axes[:, 0]:
    ax.set_ylabel("Lifetime T / T_max")
h_, l_ = axes[0, 0].get_legend_handles_labels()
fig.legend(h_, l_, loc="upper center", ncol=7, frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.95))
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(ROOT, "figures", f"fig1_fronts_same_run.{ext}"), dpi=130)
print(summ.round(4))
print(V.round(3))
print("\n".join(lines))
