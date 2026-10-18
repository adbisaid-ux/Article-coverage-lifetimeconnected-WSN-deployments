"""Revision round 1 (I10, B4): round-by-round simulation beyond first node death, and radio sensitivity.
1. Re-runs NSGA-II+SA (deterministic, same seeds) on the median-HV run of two scenarios, takes every
   distinct point of the final front, and simulates rounds with true link distances: each round the
   load-balanced minimum-hop tree is rebuilt over the alive sensors, every connected alive sensor sends
   one packet, energy is drawn per the first-order radio model, and dead sensors are removed.
   Records FND, the rounds at which connected coverage falls below 90 % and 50 % of its initial value,
   and the coverage curve.
2. Radio sensitivity: T(c)/T(1) for alternative radio parameters (the Pareto set in (CR, c) does not depend
   on them because f2 is decreasing in c for any parameters; only the mapping c -> rounds changes).
Writes results/tables_rev1/T13_round_simulation.md, T14_radio_sensitivity.md, figures/fig5_coverage_over_rounds."""
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import RADIO, Problem, bfs_hops, lifetime_from_c   # noqa: E402
from moea import nsga2_sa                                       # noqa: E402

CFG = json.load(open(os.path.join(ROOT, "experiments", "config.json")))
S = CFG["seeds"]
OUT = os.path.join(ROOT, "results", "tables_rev1")
raw = pd.read_csv(os.path.join(ROOT, "results", "raw_runs.csv"))


def simulate(xy, prob, max_rounds=20000, radio=RADIO):
    xy = np.asarray(xy, float).reshape(-1, 2)
    n = len(xy)
    E = np.full(n, radio["E0"])
    alive = np.ones(n, bool)
    k, Ee, efs = radio["k_bits"], radio["E_elec"], radio["eps_fs"]
    rs2, rc2 = prob.rs ** 2, prob.rc ** 2
    d2t = ((xy[:, None, :] - prob.targets[None, :, :]) ** 2).sum(-1) <= rs2
    curve, fnd = [], None
    for t in range(max_rounds):
        sub = np.flatnonzero(alive)
        hop_s, order_s = bfs_hops(xy[sub], prob.sink, rc2)
        conn = sub[hop_s > 0]
        cr = d2t[conn].any(0).mean() if conn.size else 0.0
        curve.append(cr)
        if conn.size == 0:
            break
        hop = np.full(n, -1)
        hop[sub] = hop_s
        order = sub[order_s]
        load = np.zeros(n)
        load[order] = 1.0
        parent = np.full(n, -1)
        for v in order[::-1]:
            if hop[v] == 1:
                continue
            cand = np.flatnonzero((hop == hop[v] - 1) & (((xy - xy[v]) ** 2).sum(1) <= rc2))
            u = cand[np.argmin(load[cand])]
            load[u] += load[v]
            parent[v] = u
        for v in order:
            p = prob.sink if parent[v] < 0 else xy[parent[v]]
            d = np.hypot(*(xy[v] - p))
            E[v] -= load[v] * k * (Ee + efs * d * d) + (load[v] - 1) * k * Ee
        dead = alive & (E <= 0)
        if dead.any() and fnd is None:
            fnd = t + 1
        alive &= ~dead
    curve = np.array(curve)
    c0 = curve[0]

    def first_below(f):
        idx = np.flatnonzero(curve < f * c0)
        return int(idx[0]) if idx.size else len(curve)
    return dict(FND=fnd, L90=first_below(0.9), L50=first_below(0.5), curve=curve)


rows, curves = [], {}
for k in ["N50_Rc12_L100", "N30_Rc8_L140"]:
    s = CFG["scenarios"][k]
    h = raw[(raw.scenario == k) & (raw.algorithm == "NSGA-II+SA")].set_index("run").HV
    r = int((h - h.median()).abs().idxmin())
    prob = Problem(s["N"], s["T"], s["L"], CFG["rs"], s["Rc"], target_seed=S["target_seed"][r])
    pop0 = np.random.default_rng(S["init_seed"][r]).uniform(0, s["L"], (CFG["NP"], prob.dim))
    res = nsga2_sa(prob, CFG["budget"], pop0, np.random.default_rng(S["algo_seed"][r]), p_sa=CFG["p_sa"])
    assert abs(prob.n_evals - CFG["budget"]) == 0
    F, X = res["F"], res["X"]
    _, idx = np.unique(np.round(F, 9), axis=0, return_index=True)
    for i in idx[np.argsort(F[idx, 0])]:
        sim = simulate(X[i], prob)
        T_model = F[i, 1] * prob.T_max
        rows.append(dict(scenario=k, run=r, CR=F[i, 0], T_model=T_model, FND_sim=sim["FND"],
                         L90_cov=sim["L90"], L50_cov=sim["L50"]))
        curves[(k, round(F[i, 0], 2), int(T_model))] = sim["curve"]
T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, "round_simulation.csv"), index=False)
open(os.path.join(OUT, "T13_round_simulation.md"), "w").write(
    "### Round-by-round simulation of every distinct point of the median NSGA-II+SA front (true distances, "
    "tree rebuilt each round over alive sensors)\n\n" + T.round(3).to_markdown(index=False) +
    "\n\nFND_sim: first node death; L90_cov / L50_cov: first round with connected coverage below 90 % / 50 % "
    "of its initial value.\n")

fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
for ax, k in zip(axes, ["N50_Rc12_L100", "N30_Rc8_L140"]):
    for (kk, cr, tm), c in curves.items():
        if kk == k:
            ax.plot(np.arange(len(c)), 100 * c, lw=1.2, label=f"CR={100 * cr:.0f}%, T={tm}")
    ax.set_xscale("symlog", linthresh=100)
    ax.set_xlabel("Round")
    ax.set_ylabel("Connected coverage (%)")
    ax.set_title(k)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(ROOT, "figures", f"fig5_coverage_over_rounds.{ext}"), dpi=130)

# radio sensitivity
sets = {"baseline (50 nJ, 10 pJ)": dict(E_elec=50e-9, eps_fs=10e-12),
        "low-power electronics (10 nJ, 10 pJ)": dict(E_elec=10e-9, eps_fs=10e-12),
        "strong amplifier (50 nJ, 100 pJ)": dict(E_elec=50e-9, eps_fs=100e-12),
        "both (10 nJ, 100 pJ)": dict(E_elec=10e-9, eps_fs=100e-12)}
rr = []
for name, p in sets.items():
    rad = dict(RADIO, **p)
    for rc in (8, 12):
        t1 = lifetime_from_c(1.0, rc, rad)
        rr.append(dict(radio=name, Rc=rc, T1=int(t1), **{f"T(c={c})/T(1)": lifetime_from_c(c, rc, rad) / t1
                                                          for c in (2, 3, 5, 10)}))
R = pd.DataFrame(rr)
open(os.path.join(OUT, "T14_radio_sensitivity.md"), "w").write(
    "### Lifetime ratio T(c)/T(1) under alternative radio parameters (E0 = 0.5 J, k = 4000 bits)\n\n" +
    R.round(3).to_markdown(index=False) + "\n")
print(T.round(3).to_string())
print(R.round(3).to_string())
