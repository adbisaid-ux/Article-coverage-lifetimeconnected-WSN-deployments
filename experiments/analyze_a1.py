"""A1 analysis: tables, tests, plausibility checks, figures (only from results/*.csv).
Usage: python3 experiments/analyze_a1.py"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG = json.load(open(os.path.join(ROOT, "experiments", "config.json")))
RES = os.path.join(ROOT, "results")
TAB = os.path.join(RES, "tables")
FIG = os.path.join(ROOT, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
MAIN = CFG["algorithms"]
PROP = "NSGA-II+SA"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#7d7d7d"]
COL = dict(zip(MAIN, PAL))
MK = dict(zip(MAIN, ["o", "s", "^", "D", "v", "P", "x"]))

raw = pd.read_csv(os.path.join(RES, "raw_runs.csv"))
keys = list(CFG["scenarios"])


def holm(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    o = np.argsort(q)
    adj = np.empty(q.size)
    run = 0.0
    for r, i in enumerate(o):
        run = max(run, (q.size - r) * q[i])
        adj[i] = min(run, 1.0)
    out[ok] = adj
    return out


def a12(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return ((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size)


def compare(df, algs, metric, ref=PROP):
    rows = []
    for k in keys:
        w = df[df.scenario == k].pivot(index="run", columns="algorithm", values=metric)[algs].dropna()
        for o in algs:
            if o == ref:
                continue
            d = w[ref] - w[o]
            p = np.nan if (d == 0).all() else stats.wilcoxon(w[ref], w[o]).pvalue
            rows.append(dict(scenario=k, comparison=f"{ref} vs {o}", metric=metric, n=len(w),
                             mean_ref=w[ref].mean(), mean_other=w[o].mean(), p=p,
                             A12=a12(w[ref], w[o]), wins=int((d > 0).sum()), losses=int((d < 0).sum())))
    st = pd.DataFrame(rows)
    st["p_holm"] = holm(st.p.values)
    st["result"] = np.where(st.p_holm < 0.05, np.where(st.mean_ref > st.mean_other, "+", "-"), "=")
    return st


def table(df, algs, metric, st, fmt, name, title, note):
    t = {}
    for k in keys:
        for a in algs:
            v = df[(df.scenario == k) & (df.algorithm == a)][metric]
            s = st[(st.scenario == k) & (st.comparison == f"{PROP} vs {a}")]
            mark = "" if a == PROP or s.empty else f" ({s.result.iloc[0]})"
            t.setdefault(k, {})[a] = f"{v.mean():{fmt}} ± {v.std(ddof=1):{fmt}}{mark}"
    tab = pd.DataFrame(t).T[algs]
    tab.index.name = "Scenario"
    tab.to_csv(os.path.join(TAB, name + ".csv"))
    open(os.path.join(TAB, name + ".md"), "w").write(f"### {title}\n\n{tab.to_markdown()}\n\n{note}\n")
    return tab


main = raw[raw.algorithm.isin(MAIN)]
st_hv = compare(main, MAIN, "HV")
st_cr = compare(main, MAIN, "max_CR")
pd.concat([st_hv, st_cr]).to_csv(os.path.join(TAB, "stats_tests.csv"), index=False)
note = (f"Mean ± sd over {CFG['runs']} runs; budget {CFG['budget']} evaluations, population {CFG['NP']}. "
        f"(+/−/=): {PROP} significantly better / worse / no difference (paired Wilcoxon signed-rank, "
        f"Holm over {len(st_hv)} comparisons, α = 0.05).")
table(main, MAIN, "HV", st_hv, ".4f", "T1_hypervolume", "Table 1. Hypervolume (objectives 1−CR, 1−T/T_max; reference point (1, 1))", note)
table(main, MAIN, "max_CR", st_cr, ".3f", "T2_max_coverage", "Table 2. Largest connected coverage on the final front", note)

# Friedman + mean ranks over all (scenario, run) blocks
w = main.pivot_table(index=["scenario", "run"], columns="algorithm", values="HV")[MAIN].dropna()
fr = stats.friedmanchisquare(*[w[a] for a in MAIN])
ranks = w.rank(axis=1, ascending=False).mean().sort_values()
with open(os.path.join(TAB, "T3_friedman_ranks.md"), "w") as fh:
    fh.write(f"### Table 3. Mean Friedman rank of HV over {len(w)} (scenario, run) blocks\n\n")
    fh.write(ranks.round(3).to_frame("mean rank").to_markdown())
    fh.write(f"\n\nFriedman χ² = {fr.statistic:.2f}, p = {fr.pvalue:.3g}.\n")

# ablation
abl_names = ["NSGA-II", "NSGA-II+SA p_sa=0.1", PROP, "NSGA-II+SA p_sa=0.5"]
abl = raw[raw.algorithm.isin(abl_names)]
if abl.algorithm.nunique() == 4:
    st_ab = compare(abl, abl_names, "HV")
    st_ab.to_csv(os.path.join(TAB, "stats_ablation.csv"), index=False)
    table(abl, abl_names, "HV", st_ab, ".4f", "T4_ablation_hv",
          "Table 4. Ablation of the sensor-agent operator (p_sa = 0, 0.1, 0.2, 0.5)",
          note.replace(f"Holm over {len(st_hv)}", f"Holm over {len(st_ab)}"))

# runtime
rt = main.groupby("algorithm").time_s.agg(["mean", "std"]).loc[MAIN].round(2)
rt.to_csv(os.path.join(TAB, "T5_runtime.csv"))
open(os.path.join(TAB, "T5_runtime.md"), "w").write("### Table 5. Wall time per run (s), all scenarios\n\n" + rt.to_markdown() + "\n")

# validation of the lifetime model
val = pd.read_csv(os.path.join(RES, "validation.csv"))
val["ratio_exact"] = val.T_tree / val.T_exact_flow
val["ratio_true"] = val.FND_true_dist / np.floor(val.T_tree)
vs = val.groupby("scenario")[["ratio_exact", "ratio_true"]].describe(percentiles=[.05, .5, .95])
vs.to_csv(os.path.join(TAB, "T6_model_validation.csv"))
rho = stats.spearmanr(val.T_tree, val.FND_true_dist)
with open(os.path.join(TAB, "T6_model_validation.md"), "w") as fh:
    fh.write("### Table 6. Lifetime model validation on final-front points\n\n")
    g = val.groupby("scenario").agg(n=("ratio_exact", "size"),
                                    tree_over_flow_median=("ratio_exact", "median"),
                                    tree_over_flow_min=("ratio_exact", "min"),
                                    fnd_true_over_tree_median=("ratio_true", "median"),
                                    fnd_true_over_tree_min=("ratio_true", "min"))
    fh.write(g.round(3).to_markdown())
    fh.write(f"\n\nT_tree / T_flow ≤ 1 by construction (load-balanced tree vs optimal flow routing). "
             f"FND with true distances / floor(T_tree) ≥ 1 expected (all links ≤ Rc). "
             f"Spearman(T_tree, FND_true) = {rho.statistic:.4f} over {len(val)} points.\n")

# ---------------------------------------------------------------- plausibility (rule 6.4)
rep = ["# A1 : contrôles de vraisemblance automatiques\n"]
exp_n = len(keys) * CFG["runs"] * (len(MAIN) + len(CFG["ablation"]))
rep.append(f"- Exécutions : {len(raw)} / {exp_n}.")
over = raw.n_evals - raw.budget
rep.append(f"- Évaluations au-delà du budget : max {over.max()} ({100 * over.max() / CFG['budget']:.2f} %), "
           f"par algorithme : {over.groupby(raw.algorithm).max().to_dict()}.")
rep.append(f"- Évaluations en deçà du budget : {int((over < 0).sum())} exécution(s).")
rep.append(f"- HV hors [0, 1] : {int(((raw.HV < 0) | (raw.HV > 1)).sum())} ; CR hors [0, 1] : "
           f"{int(((raw.max_CR < 0) | (raw.max_CR > 1 + 1e-9)).sum())}.")
rep.append(f"- Validation : T_tree > T_flow au-delà de la tolérance de bisection (1e-3, impossible) : "
           f"{int((val.ratio_exact > 1 + 1e-3).sum())} (dépassement max {val.ratio_exact.max() - 1:.1e}, "
           f"dû à c* calculé par excès à 1e-3 près) ; points où l'arbre perd plus de 0,1 % : "
           f"{100 * (val.ratio_exact < 0.999).mean():.1f} % ; "
           f"FND (vraies distances) < floor(T_tree) : {int((val.ratio_true < 1 - 1e-9).sum())}.")
g = raw.groupby(["scenario", "algorithm"]).HV.std()
z = g[g < 1e-9]
rep.append(f"- Écart-type de HV nul : {len(z)} cellule(s) {list(z.index)}.")
g = raw.groupby(["scenario", "algorithm"]).max_Tn.agg(["min", "max"])
rep.append(f"- Cellules où T/T_max maximal vaut 1 dans toutes les exécutions : {int((g['min'] > 1 - 1e-9).sum())} "
           f"sur {len(g)} (attendu : un nœud par saut au puits suffit à atteindre c = 1).")
open(os.path.join(RES, "plausibility_report.md"), "w").write("\n".join(rep) + "\n")

# ---------------------------------------------------------------- figures
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#dddddd", "grid.linewidth": 0.6})
fronts = pd.read_csv(os.path.join(RES, "fronts.csv"))
fig, axes = plt.subplots(2, 4, figsize=(13, 6), sharey=True)
for ax, k in zip(axes.ravel(), keys):
    for a in MAIN:
        d = main[(main.scenario == k) & (main.algorithm == a)]
        r = int(d.iloc[(d.HV - d.HV.median()).abs().argmin()].run)    # median-HV run
        f = fronts[(fronts.scenario == k) & (fronts.algorithm == a) & (fronts.run == r)].sort_values("CR")
        ax.plot(f.CR * 100, f.Tn, color=COL[a], marker=MK[a], ms=3, lw=1.2, label=a)
    ax.set_title(k)
    ax.set_xlabel("Connected coverage (%)")
for ax in axes[:, 0]:
    ax.set_ylabel("Lifetime T / T_max")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=len(MAIN), frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for e in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig1_fronts_median_run.{e}"), dpi=200)

hist = pd.read_csv(os.path.join(RES, "hv_history.csv"))
fig, axes = plt.subplots(2, 4, figsize=(13, 6), sharex=True)
grid = np.linspace(CFG["budget"] / 20, CFG["budget"], 20)
for ax, k in zip(axes.ravel(), keys):
    for a in MAIN:
        H = hist[(hist.scenario == k) & (hist.algorithm == a)]
        M = []
        for r, g in H.groupby("run"):
            g = g.sort_values("fe")
            idx = np.searchsorted(g.fe.values, grid, side="right") - 1
            M.append(g.HV.values[np.maximum(idx, 0)])
        ax.plot(grid, np.median(M, 0), color=COL[a], lw=1.5, label=a)
    ax.set_title(k)
for ax in axes[1]:
    ax.set_xlabel("Evaluations")
for ax in axes[:, 0]:
    ax.set_ylabel("HV (median of 30 runs)")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=len(MAIN), frameon=False)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for e in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig2_hv_convergence.{e}"), dpi=200)

if abl.algorithm.nunique() == 4:
    fig, axes = plt.subplots(1, 8, figsize=(14, 3), sharey=False)
    for ax, k in zip(axes, keys):
        data = [abl[(abl.scenario == k) & (abl.algorithm == a)].HV.values for a in abl_names]
        ax.boxplot(data, widths=0.6, medianprops=dict(color="black"))
        ax.set_xticks(range(1, 5), ["0", "0.1", "0.2", "0.5"])
        ax.set_title(k, fontsize=7)
        ax.set_xlabel("p_sa")
    axes[0].set_ylabel("HV")
    fig.tight_layout()
    for e in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, f"fig3_ablation.{e}"), dpi=200)

fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
ax[0].scatter(val.T_exact_flow, val.T_tree, s=4, alpha=0.4, color=PAL[0])
m = max(val.T_exact_flow.max(), val.T_tree.max())
ax[0].plot([0, m], [0, m], color="#555", lw=1)
ax[0].set_xlabel("Lifetime, optimal flow routing (rounds)")
ax[0].set_ylabel("Lifetime, balanced tree (rounds)")
ax[1].scatter(val.T_tree, val.FND_true_dist, s=4, alpha=0.4, color=PAL[2])
m = max(val.T_tree.max(), val.FND_true_dist.max())
ax[1].plot([0, m], [0, m], color="#555", lw=1)
ax[1].set_xlabel("Model lifetime T_tree (rounds)")
ax[1].set_ylabel("Simulated FND, true distances")
fig.tight_layout()
for e in ("png", "pdf"):
    fig.savefig(os.path.join(FIG, f"fig4_model_validation.{e}"), dpi=200)
print("done")
