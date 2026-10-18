"""Aggregation-regime campaign (results_agg): statistics, controls, trade-off, cross-regime evaluation.
Paired two-sided Wilcoxon signed-rank tests with Holm correction within each stated family.
Writes results_agg/tables/*.md|csv."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import RADIO, lifetime_agg, lifetime_from_c   # noqa: E402
from moea import hv2d                                    # noqa: E402

RES = os.path.join(ROOT, "results_agg")
OUT = os.path.join(RES, "tables")
os.makedirs(OUT, exist_ok=True)
raw = pd.read_csv(os.path.join(RES, "raw_runs.csv"))
fr = pd.read_csv(os.path.join(RES, "fronts.csv"))
val = pd.read_csv(os.path.join(RES, "validation.csv"))
main = pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv"))
boucle = pd.read_csv(os.path.join(ROOT, "results_boucle", "raw_runs.csv"))
AG = raw[raw.scenario.str.startswith("A_")].copy()
XR = raw[~raw.scenario.str.startswith("A_")].copy()
AG["sc"] = AG.scenario.str[2:]
SC = ["N30_Rc8_L100", "N30_Rc8_L140", "N30_Rc12_L100", "N30_Rc12_L140",
      "N50_Rc8_L100", "N50_Rc8_L140", "N50_Rc12_L100", "N50_Rc12_L140"]
lines = []


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(run, 1.0)
    return adj


def a12(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return ((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (len(x) * len(y))


def compare(df, base, others, metric, family):
    rows = []
    for k in SC:
        w = df[df.sc == k].pivot(index="run", columns="algorithm", values=metric)
        for o in others:
            rows.append(dict(scenario=k, metric=metric, a=base, b=o, ma=w[base].mean(), mb=w[o].mean(),
                             p=stats.wilcoxon(w[base], w[o]).pvalue if (w[base] != w[o]).any() else 1.0,
                             A12=a12(w[base], w[o])))
    r = pd.DataFrame(rows)
    r["p_holm"] = holm(r.p)
    r["res"] = np.where(r.p_holm < 0.05, np.where(r.ma > r.mb, "+", "-"), "=")
    r["gain%"] = 100 * (r.ma / r.mb - 1)
    r["family"] = family
    return r


# 0. reproducibility: XR replicates the main campaign (NSGA-II, NSGA-II+SA) and E1 (NSGA-II+REPAIR)
ref = pd.concat([main, boucle[boucle.algorithm == "NSGA-II+REPAIR"]])
chk = XR.merge(ref, on=["scenario", "algorithm", "run"], suffixes=("", "_ref"))
lines.append(f"Reproducibility: {len(chk)} XR runs matched to earlier runs; max |HV - HV_ref| = "
             f"{(chk.HV - chk.HV_ref).abs().max():.3e}")

# 1. HV and largest coverage, NSGA-II+SA against every other algorithm (family: 48 tests per metric)
algs = sorted(AG.algorithm.unique())
others = [a for a in ["NSGA-II", "NSGA-II+REPAIR", "NSGA-II+RND", "NSGA-III", "SMS-EMOA", "AGE-MOEA-II"] if a in algs]
S1 = pd.concat([compare(AG, "NSGA-II+SA", others, m, "SA vs all, per metric") for m in ["HV", "max_CR"]])
S1.to_csv(os.path.join(OUT, "stats_agg.csv"), index=False)
g = AG.groupby(["sc", "algorithm"]).HV.agg(["mean", "std"]).unstack()
tab = g["mean"].round(4).astype(str) + " ± " + g["std"].round(4).astype(str)
lines.append("### Aggregation regime: HV mean ± sd over 30 runs\n\n" + tab.loc[SC].to_markdown())
c = AG.groupby(["sc", "algorithm"]).max_CR.mean().unstack()
lines.append("### Aggregation regime: largest connected coverage, mean over 30 runs\n\n" + c.loc[SC].round(3).to_markdown())
for m in ["HV", "max_CR"]:
    q = S1[S1.metric == m]
    lines.append(f"{m}: NSGA-II+SA vs others, Holm over {len(q)}: " +
                 ", ".join(f"{o}: {(q[q.b == o].res == '+').sum()}+/{(q[q.b == o].res == '=').sum()}=/"
                           f"{(q[q.b == o].res == '-').sum()}- (gain {q[q.b == o]['gain%'].min():.1f} to "
                           f"{q[q.b == o]['gain%'].max():.1f} %)" for o in others))

# 2. controls: pairwise among SA, REPAIR, RND (family: 24 tests per metric)
rows = []
for m in ["HV", "max_CR"]:
    rr = []
    for a, b in [("NSGA-II+SA", "NSGA-II+REPAIR"), ("NSGA-II+SA", "NSGA-II+RND"), ("NSGA-II+REPAIR", "NSGA-II+RND")]:
        rr.append(compare(AG, a, [b], m, "controls"))
    rr = pd.concat(rr)
    rr["p_holm"] = holm(rr.p)
    rr["res"] = np.where(rr.p_holm < 0.05, np.where(rr.ma > rr.mb, "+", "-"), "=")
    rows.append(rr)
C = pd.concat(rows)
C.to_csv(os.path.join(OUT, "stats_agg_controls.csv"), index=False)
lines.append("### Controls in the aggregation regime (Holm over 24 tests per metric)\n\n" +
             C[["scenario", "metric", "a", "b", "ma", "mb", "p_holm", "res", "A12", "gain%"]].round(4).to_markdown(index=False))

# 3. trade-off: median largest coverage per in-degree level d (aggregation) for NSGA-II+SA
rc_of = {k: (8.0 if "Rc8" in k else 12.0) for k in SC}
f = fr[fr.scenario.str.startswith("A_") & (fr.algorithm == "NSGA-II+SA")].copy()
f["sc"] = f.scenario.str[2:]
rows = []
for k in SC:
    rc = rc_of[k]
    ff = f[f.sc == k]
    for d in range(0, 6):
        tn = lifetime_agg(d, rc) / lifetime_agg(0, rc)
        best = ff[ff.Tn >= tn - 1e-9].groupby("run").CR.max().reindex(range(30)).fillna(0)
        rows.append(dict(scenario=k, d=d, Tn=tn, median_CR=100 * best.median()))
TD = pd.DataFrame(rows)
lines.append("### Aggregation regime, NSGA-II+SA: median largest coverage (%) with in-degree at most d\n\n" +
             TD.pivot(index="scenario", columns="d", values="median_CR").loc[SC].round(0).to_markdown())
lines.append("Lifetime ratio T_agg(d)/T_agg(0): " + ", ".join(
    f"Rc={rc:g}: " + " / ".join(f"{lifetime_agg(d, rc) / lifetime_agg(0, rc):.3f}" for d in range(6)) for rc in (8.0, 12.0)))

# 4. cross-regime: what a front optimised for one regime is worth under the other
#    For run r of scenario k: HV under regime X of (i) the front optimised for X and (ii) the front optimised
#    for the other regime, re-evaluated under X and filtered to its non-dominated points.
rows = []
for k in SC:
    for alg in ["NSGA-II+SA", "NSGA-II"]:
        for r in range(30):
            a = fr[(fr.scenario == "A_" + k) & (fr.algorithm == alg) & (fr.run == r)]
            b = fr[(fr.scenario == k) & (fr.algorithm == alg) & (fr.run == r)]
            if a.empty or b.empty:
                continue
            hv_agg_own = hv2d(1 - a[["CR", "Tn"]].values)
            hv_agg_from_raw = hv2d(1 - b[["CR", "Tn_other"]].values)
            hv_raw_own = hv2d(1 - b[["CR", "Tn"]].values)
            hv_raw_from_agg = hv2d(1 - a[["CR", "Tn_other"]].values)
            rows.append(dict(scenario=k, algorithm=alg, run=r, hv_agg_own=hv_agg_own, hv_agg_from_raw=hv_agg_from_raw,
                             hv_raw_own=hv_raw_own, hv_raw_from_agg=hv_raw_from_agg))
X = pd.DataFrame(rows)
X.to_csv(os.path.join(OUT, "cross_regime.csv"), index=False)
rows = []
for (k, alg), d in X.groupby(["scenario", "algorithm"]):
    rows.append(dict(scenario=k, algorithm=alg,
                     agg_loss_pct=100 * (1 - d.hv_agg_from_raw.mean() / d.hv_agg_own.mean()),
                     p_agg=stats.wilcoxon(d.hv_agg_own, d.hv_agg_from_raw).pvalue,
                     raw_loss_pct=100 * (1 - d.hv_raw_from_agg.mean() / d.hv_raw_own.mean()),
                     p_raw=stats.wilcoxon(d.hv_raw_own, d.hv_raw_from_agg).pvalue))
XS = pd.DataFrame(rows)
XS["p_agg_holm"] = holm(XS.p_agg)
XS["p_raw_holm"] = holm(XS.p_raw)
XS.to_csv(os.path.join(OUT, "cross_regime_stats.csv"), index=False)
lines.append("### Cross-regime HV loss (% of the HV of the front optimised for the evaluating regime); "
             "Holm over 16 tests per column\n\n" + XS.round(4).to_markdown(index=False))

# 5. validation of the aggregation lifetime
v = val[val.scenario.str.startswith("A_")]
lines.append(f"Validation (aggregation): {len(v)} front points; T_model == T_exact (max-flow) in "
             f"{(np.abs(v.T_tree - v.T_exact_flow) < 1e-6).mean() * 100:.1f} %; true-distance FND >= model in "
             f"{(v.FND_true_dist >= v.T_tree - 1e-6).mean() * 100:.1f} %; median ratio FND_true/T_model "
             f"{(v.FND_true_dist / v.T_tree).median():.3f}")
lines.append("Run time (s), mean per algorithm: " + AG.groupby("algorithm").time_s.mean().round(1).to_dict().__repr__())
lines.append(f"Budget check: n_evals == budget in {(raw.n_evals == raw.budget).mean() * 100:.0f} % of runs")
open(os.path.join(OUT, "summary_agg.md"), "w").write("\n\n".join(lines) + "\n")
print("\n\n".join(lines))

# 6. figure: fronts of NSGA-II+SA under the two regimes (run with median HV under aggregation), and
#    largest coverage reachable at lifetime >= T_max/3 in each regime
import matplotlib   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

fig, axes = plt.subplots(2, 4, figsize=(14, 6.2), sharey=True)
for ax, k in zip(axes.ravel(), SC):
    h = AG[(AG.sc == k) & (AG.algorithm == "NSGA-II+SA")].set_index("run").HV
    r = int((h - h.median()).abs().idxmin())
    for sc, lab, col, mk in [("A_" + k, "optimised with aggregation", "tab:blue", "o"),
                             (k, "optimised without aggregation", "tab:red", "s")]:
        p = fr[(fr.scenario == sc) & (fr.algorithm == "NSGA-II+SA") & (fr.run == r)]
        p = p.drop_duplicates(["CR", "Tn"]).sort_values("CR")
        ax.step(100 * p.CR, p.Tn, where="pre", color=col, lw=1.2)
        ax.plot(100 * p.CR, p.Tn, mk, color=col, ms=3.5, label=lab)
    ax.set_title(k, fontsize=9)
    ax.set_xlabel("Connected coverage (%)", fontsize=8)
    ax.grid(alpha=0.3)
axes[0, 0].set_ylabel("Normalised lifetime $f_2$ (own regime)")
axes[1, 0].set_ylabel("Normalised lifetime $f_2$ (own regime)")
axes[0, 0].legend(fontsize=7)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(ROOT, "figures", f"fig6_fronts_two_regimes.{ext}"), dpi=130)

# coverage at lifetime >= T_max/3 (c <= 2 without aggregation, d <= 2 with aggregation), NSGA-II+SA, median over runs
rows = []
for k in SC:
    rc = rc_of[k]
    ra = fr[(fr.scenario == k) & (fr.algorithm == "NSGA-II+SA")]
    ag = fr[(fr.scenario == "A_" + k) & (fr.algorithm == "NSGA-II+SA")]
    t_raw = lifetime_from_c(2.0, rc) / lifetime_from_c(1.0, rc)
    t_agg = lifetime_agg(2, rc) / lifetime_agg(0, rc)
    cr_raw = ra[ra.Tn >= t_raw - 1e-9].groupby("run").CR.max().reindex(range(30)).fillna(0)
    cr_agg = ag[ag.Tn >= t_agg - 1e-9].groupby("run").CR.max().reindex(range(30)).fillna(0)
    rows.append(dict(scenario=k, ratio_raw_c2=t_raw, ratio_agg_d2=t_agg, median_CR_raw=100 * cr_raw.median(),
                     median_CR_agg=100 * cr_agg.median()))
Q = pd.DataFrame(rows)
Q.to_csv(os.path.join(OUT, "coverage_at_third_lifetime.csv"), index=False)
open(os.path.join(OUT, "summary_agg.md"), "a").write(
    "\n### Largest coverage (median %) with lifetime >= about T_max/3 (c <= 2 raw, d <= 2 aggregation), NSGA-II+SA\n\n"
    + Q.round(3).to_markdown(index=False) + "\n")
print(Q.round(3).to_string())
