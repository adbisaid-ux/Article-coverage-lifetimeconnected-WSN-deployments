"""Sensitivity to a per-round sensing energy E_s (review round 4, item 10).

Every sensor spends E_s per round in addition to radio energy, so
    T(c) = E0 / (c e_tx + (c - 1) e_rx + E_s)      without aggregation,
    T_agg(d) = E0 / (e_tx + d e_rx + E_s)          with aggregation.
Both are decreasing in c (resp. d), so f2 is a monotone transform of the E_s = 0 objective: the Pareto set in
decision space is unchanged and the stored final fronts can be re-scored exactly. Tn -> c (or d) is inverted
from the stored normalised lifetime, then f2 is recomputed with E_s and the hypervolume of every run is recomputed.
Writes results_sensing/: ratios.csv, hv_runs.csv, stats_hv.csv, coverage_equal_lifetime.csv, summary.md."""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import radio_costs   # noqa: E402
from moea import hv2d           # noqa: E402

OUT = os.path.join(ROOT, "results_sensing")
os.makedirs(OUT, exist_ok=True)
LEVELS = [0.0, 0.5, 1.0, 2.0, 5.0]          # E_s in units of e_rx = k E_elec (receiving one packet)


def rc_of(scenario):
    return float(scenario.split("_Rc")[1].split("_")[0])


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for i, j in enumerate(o):
        run = max(run, (len(p) - i) * p[j])
        adj[j] = min(run, 1.0)
    return adj


def rescore(Tn, rc, regime, es_units):
    """Normalised lifetime under sensing energy E_s = es_units * e_rx, from the stored E_s = 0 value."""
    e_tx, e_rx = radio_costs(rc)
    Es = es_units * e_rx
    Tn = np.asarray(Tn, float)
    out = np.zeros_like(Tn)
    ok = Tn > 0
    if regime == "raw":
        c = (e_tx / Tn[ok] + e_rx) / (e_tx + e_rx)
        out[ok] = (e_tx + Es) / (c * e_tx + (c - 1) * e_rx + Es)
    else:
        d = (e_tx / Tn[ok] - e_tx) / e_rx
        out[ok] = (e_tx + Es) / (e_tx + d * e_rx + Es)
    return out


# 1. analytic ratios
rows = []
for rc in (8.0, 12.0):
    e_tx, e_rx = radio_costs(rc)
    for u in LEVELS:
        Es = u * e_rx
        T = lambda c: 1.0 / (c * e_tx + (c - 1) * e_rx + Es)          # noqa: E731
        Ta = lambda d: 1.0 / (e_tx + d * e_rx + Es)                    # noqa: E731
        rows.append(dict(Rc=rc, Es_over_erx=u, T2_over_T1=T(2) / T(1), T6_over_T1=T(6) / T(1),
                         Tagg2_over_Tagg0=Ta(2) / Ta(0),
                         max_d_with_Tagg_ge_T2=int(np.floor(1 + e_tx / e_rx + 1e-12))))
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "ratios.csv"), index=False)

# 2. re-scored hypervolumes of every run (main campaign, raw regime; aggregation campaign)
sources = [("results", "raw"), ("results_agg", "agg")]
hv_rows = []
for folder, default_regime in sources:
    fr = pd.read_csv(os.path.join(ROOT, folder, "fronts.csv"))
    for (sc, alg, run), g in fr.groupby(["scenario", "algorithm", "run"]):
        regime = "agg" if sc.startswith("A_") else "raw"
        rc = rc_of(sc)
        for u in LEVELS:
            F = np.c_[g.CR.values, rescore(g.Tn.values, rc, regime, u)]
            hv_rows.append(dict(source=folder, scenario=sc, regime=regime, algorithm=alg, run=run,
                                Es_over_erx=u, HV=hv2d(1.0 - F)))
H = pd.DataFrame(hv_rows)
H.to_csv(os.path.join(OUT, "hv_runs.csv"), index=False)

# sanity check: E_s = 0 reproduces the stored hypervolumes
check = []
for folder, _ in sources:
    raw = pd.read_csv(os.path.join(ROOT, folder, "raw_runs.csv"))
    m = H[(H.source == folder) & (H.Es_over_erx == 0)].merge(raw, on=["scenario", "algorithm", "run"],
                                                             suffixes=("", "_stored"))
    check.append(np.abs(m.HV - m.HV_stored).max())
print("max |HV(E_s=0) - stored HV|:", check)

# 3. NSGA-II+SA vs each competitor at every E_s level: paired Wilcoxon, Holm within (source, level)
st = []
for (src, u), d in H.groupby(["source", "Es_over_erx"]):
    rows = []
    for sc, ds in d.groupby("scenario"):
        w = ds.pivot(index="run", columns="algorithm", values="HV")
        for o in w.columns:
            if o == "NSGA-II+SA" or o.startswith("NSGA-II+SA "):
                continue
            p = stats.wilcoxon(w["NSGA-II+SA"], w[o]).pvalue if (w["NSGA-II+SA"] != w[o]).any() else 1.0
            rows.append(dict(source=src, Es_over_erx=u, scenario=sc, other=o, sa=w["NSGA-II+SA"].mean(),
                             mo=w[o].mean(), p=p))
    r = pd.DataFrame(rows)
    r["p_holm"] = holm(r.p)
    r["res"] = np.where(r.p_holm < 0.05, np.where(r.sa > r.mo, "+", "-"), "=")
    st.append(r)
S = pd.concat(st)
S.to_csv(os.path.join(OUT, "stats_hv.csv"), index=False)

# 4. coverage at equal lifetime across regimes (lifetime of a c = 2 tree without aggregation)
fr_raw = pd.read_csv(os.path.join(ROOT, "results_agg", "fronts.csv"))
cov = []
for sc in sorted(s for s in fr_raw.scenario.unique() if not s.startswith("A_")):
    rc = rc_of(sc)
    e_tx, e_rx = radio_costs(rc)
    ra = fr_raw[(fr_raw.scenario == sc) & (fr_raw.algorithm == "NSGA-II+SA")]
    ag = fr_raw[(fr_raw.scenario == "A_" + sc) & (fr_raw.algorithm == "NSGA-II+SA")]
    for u in LEVELS:
        Es = u * e_rx
        level = (e_tx + Es) / (2 * e_tx + e_rx + Es)                   # T(2)/T(1) with sensing energy
        tr = rescore(ra.Tn.values, rc, "raw", u)
        ta = rescore(ag.Tn.values, rc, "agg", u)
        a = ra.assign(t=tr)[tr >= level - 1e-9].groupby("run").CR.max().reindex(range(30)).fillna(0)
        b = ag.assign(t=ta)[ta >= level - 1e-9].groupby("run").CR.max().reindex(range(30)).fillna(0)
        cov.append(dict(scenario=sc, Es_over_erx=u, lifetime_fraction=level, median_CR_raw=100 * a.median(),
                        median_CR_agg=100 * b.median(), ratio=b.median() / a.median() if a.median() > 0 else np.nan))
C = pd.DataFrame(cov)
C.to_csv(os.path.join(OUT, "coverage_equal_lifetime.csv"), index=False)

with open(os.path.join(OUT, "summary.md"), "w") as f:
    f.write("# Sensitivity to sensing energy\n\n## Analytic ratios\n\n" + R.round(3).to_markdown(index=False))
    f.write("\n\n## NSGA-II+SA vs competitors (HV, Holm within source and level): counts of +/=/-\n\n")
    f.write(S.groupby(["source", "Es_over_erx", "res"]).size().unstack(fill_value=0).to_markdown())
    f.write("\n\n## Coverage at the lifetime of a c = 2 tree (median %, NSGA-II+SA)\n\n" + C.round(3).to_markdown(index=False))
    f.write(f"\n\nSanity check, max |HV(E_s=0) - stored HV| per source: {check}\n")
print(open(os.path.join(OUT, "summary.md")).read())
