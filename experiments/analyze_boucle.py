"""/boucle analysis: E1 (control operators), E2 (SA-op in other MOEAs), E3 (robustness scenarios).
Reads results/raw_runs.csv (main campaign) and results_boucle/raw_runs.csv.
Writes results_boucle/tables/*.md|csv, results_boucle/plausibility_report.md and figures_boucle/."""
import json
import os

import matplotlib
import numpy as np
import pandas as pd
from scipy import stats

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG = json.load(open(os.path.join(ROOT, "experiments", "config_boucle.json")))
OUT = os.path.join(ROOT, "results_boucle")
TAB = os.path.join(OUT, "tables")
FIG = os.path.join(ROOT, "figures_boucle")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
main = pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv"))
new = pd.read_csv(os.path.join(OUT, "raw_runs.csv"))
raw = pd.concat([main, new], ignore_index=True)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for r, i in enumerate(o):
        run = max(run, (p.size - r) * p[i])
        adj[i] = min(run, 1.0)
    return adj


def a12(x, y):
    x, y = np.asarray(x), np.asarray(y)
    return ((x[:, None] > y[None, :]).sum() + 0.5 * (x[:, None] == y[None, :]).sum()) / (x.size * y.size)


def compare(pairs, keys, metric):
    """pairs: list of (a, b); paired Wilcoxon a vs b per scenario, Holm over all."""
    rows = []
    for k in keys:
        d = raw[raw.scenario == k]
        for a, b in pairs:
            w = d[d.algorithm.isin([a, b])].pivot(index="run", columns="algorithm", values=metric)[[a, b]].dropna()
            diff = w[a] - w[b]
            p = 1.0 if (diff == 0).all() else stats.wilcoxon(w[a], w[b]).pvalue
            rows.append(dict(scenario=k, a=a, b=b, metric=metric, n=len(w), mean_a=w[a].mean(),
                             sd_a=w[a].std(ddof=1), mean_b=w[b].mean(), sd_b=w[b].std(ddof=1),
                             gain_pct=100 * (w[a].mean() / w[b].mean() - 1), p=p, A12=a12(w[a], w[b])))
    st = pd.DataFrame(rows)
    st["p_holm"] = holm(st.p.values)
    st["result"] = np.where(st.p_holm < 0.05, np.where(st.mean_a > st.mean_b, "+", "-"), "=")
    return st


def cell(v):
    return f"{v.mean():.4f} ± {v.std(ddof=1):.4f}"


rep = ["# A1 /boucle : contrôles de vraisemblance automatiques", ""]
exp_runs = {"E12": len(CFG["scenarios_E12"]) * len(CFG["algorithms_E12"]) * CFG["runs"],
            "E3": len(CFG["scenarios_E3"]) * len(CFG["algorithms_E3"]) * CFG["runs"]}
rep.append(f"- Exécutions : {len(new)} / {sum(exp_runs.values())} attendues.")
over = new.n_evals - new.budget
rep.append(f"- Évaluations hors budget : max {over.max()}, min {over.min()}.")
rep.append(f"- HV hors [0, 1] : {int(((new.HV < 0) | (new.HV > 1)).sum())} ; "
           f"CR hors [0, 1] : {int(((new.max_CR < 0) | (new.max_CR > 1)).sum())}.")
g = new.groupby(["scenario", "algorithm"]).HV.std()
rep.append(f"- Écart-type de HV nul : {int((g < 1e-12).sum())} cellule(s).")

K8 = CFG["scenarios_E12"]
# ---------------------------------------------------------------- E1: control operators
e1_algs = ["NSGA-II", "NSGA-II+REPAIR", "NSGA-II+RND", "NSGA-II+SA"]
st1 = compare([("NSGA-II+SA", "NSGA-II+RND"), ("NSGA-II+SA", "NSGA-II+REPAIR"),
               ("NSGA-II+RND", "NSGA-II"), ("NSGA-II+REPAIR", "NSGA-II")], K8, "HV")
st1.to_csv(os.path.join(TAB, "stats_E1.csv"), index=False)
t = {}
for k in K8:
    for al in e1_algs:
        v = raw[(raw.scenario == k) & (raw.algorithm == al)].HV
        t.setdefault(k, {})[al] = cell(v)
    for a, b in [("NSGA-II+SA", "NSGA-II+RND"), ("NSGA-II+SA", "NSGA-II+REPAIR")]:
        s = st1[(st1.scenario == k) & (st1.a == a) & (st1.b == b)].iloc[0]
        t[k][f"SA vs {b[8:]}"] = f"{s.result} (A12 {s.A12:.2f})"
tab = pd.DataFrame(t).T
tab.index.name = "Scenario"
tab.to_csv(os.path.join(TAB, "E1_control_operators.csv"))
open(os.path.join(TAB, "E1_control_operators.md"), "w").write(
    "### E1. Control operators in NSGA-II (p = 0.2): HV mean ± sd over 30 runs\n\n" + tab.to_markdown() +
    "\n\nREPAIR: projection onto the parent disc and attachment of disconnected sensors, no move. "
    "RND: same, redundant sensors move by r_s in a random direction. Paired Wilcoxon, Holm over "
    f"{len(st1)} comparisons (incl. RND vs NSGA-II and REPAIR vs NSGA-II, see stats_E1.csv).\n")
for (a, b), s in st1.groupby(["a", "b"]):
    rep.append(f"- E1 {a} vs {b} : {''.join(s.result)} ; gain HV {s.gain_pct.min():.1f} à {s.gain_pct.max():.1f} %.")

# ---------------------------------------------------------------- E2: SA-op in other MOEAs
pairs2 = [("NSGA-III+SA", "NSGA-III"), ("SMS-EMOA+SA", "SMS-EMOA")]
st2 = pd.concat([compare(pairs2, K8, "HV"), compare(pairs2, K8, "max_CR")])
st2.to_csv(os.path.join(TAB, "stats_E2.csv"), index=False)
t = {}
for k in K8:
    for a, b in pairs2:
        for m in ["HV", "max_CR"]:
            s = st2[(st2.scenario == k) & (st2.a == a) & (st2.metric == m)].iloc[0]
            t.setdefault(k, {})[f"{b} {m}"] = f"{s.mean_b:.4f} → {s.mean_a:.4f} ({s.result})"
tab = pd.DataFrame(t).T
tab.index.name = "Scenario"
tab.to_csv(os.path.join(TAB, "E2_saop_in_other_moeas.csv"))
open(os.path.join(TAB, "E2_saop_in_other_moeas.md"), "w").write(
    "### E2. SA-op (p = 0.2, after polynomial mutation) inside pymoo NSGA-III and SMS-EMOA: "
    "mean without → with SA-op\n\n" + tab.to_markdown() +
    "\n\n(+/−/=): with SA-op significantly better / worse / no difference (paired Wilcoxon, Holm over "
    "the 16 comparisons of each metric).\n")
for (a, m), s in st2.groupby(["a", "metric"]):
    rep.append(f"- E2 {a} ({m}) : {''.join(s.result)} ; gain {s.gain_pct.min():.1f} à {s.gain_pct.max():.1f} %.")

# ---------------------------------------------------------------- E3: robustness scenarios
K3 = list(CFG["scenarios_E3"])
algs = CFG["algorithms_E3"]
others = [a for a in algs if a != "NSGA-II+SA"]
for m, fmt, name in [("HV", ".4f", "E3_hypervolume"), ("max_CR", ".3f", "E3_max_coverage")]:
    st3 = compare([("NSGA-II+SA", o) for o in others], K3, m)
    st3.to_csv(os.path.join(TAB, f"stats_{name}.csv"), index=False)
    t = {}
    for k in K3:
        for al in algs:
            v = raw[(raw.scenario == k) & (raw.algorithm == al)][m]
            s = st3[(st3.scenario == k) & (st3.b == al)]
            mark = "" if s.empty else f" ({s.result.iloc[0]})"
            t.setdefault(k, {})[al] = f"{v.mean():{fmt}} ± {v.std(ddof=1):{fmt}}{mark}"
    tab = pd.DataFrame(t).T[algs]
    tab.index.name = "Scenario"
    tab.to_csv(os.path.join(TAB, name + ".csv"))
    open(os.path.join(TAB, name + ".md"), "w").write(
        f"### E3. {m} on the robustness scenarios, mean ± sd over 30 runs\n\n" + tab.to_markdown() +
        f"\n\n(+/−/=): NSGA-II+SA significantly better / worse / no difference (paired Wilcoxon, Holm over "
        f"{len(st3)} comparisons).\n")
    rep.append(f"- E3 {m} : NSGA-II+SA vs autres : {''.join(st3.result)} ; A12 min {st3.A12.min():.2f}.")
    best = st3.loc[st3.groupby("scenario").mean_b.idxmax()]
    for _, s in best.iterrows():
        rep.append(f"  - {s.scenario} : meilleur concurrent {s.b}, gain {s.gain_pct:.1f} %.")
rt = new.groupby("algorithm").time_s.agg(["mean", "std"]).round(2)
open(os.path.join(TAB, "runtime_boucle.md"), "w").write(rt.to_markdown())

# ---------------------------------------------------------------- figure: E1 + E2 boxplots
fig, axes = plt.subplots(1, len(K8), figsize=(20, 3.6), sharey=False)
labs = ["NSGA-II", "REPAIR", "RND", "SA"]
for ax, k in zip(axes, K8):
    d = raw[raw.scenario == k]
    ax.boxplot([d[d.algorithm == a].HV for a in e1_algs], tick_labels=labs, widths=0.6)
    ax.set_title(k, fontsize=9)
    ax.tick_params(axis="x", labelrotation=45, labelsize=8)
    ax.grid(alpha=0.3)
axes[0].set_ylabel("HV")
fig.tight_layout()
fig.savefig(os.path.join(FIG, "figB1_control_operators.pdf"))
fig.savefig(os.path.join(FIG, "figB1_control_operators.png"), dpi=130)
plt.close(fig)

open(os.path.join(OUT, "plausibility_report.md"), "w").write("\n".join(rep) + "\n")
print("\n".join(rep))
