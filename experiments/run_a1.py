"""A1 campaign. Usage: python3 experiments/run_a1.py [--workers 4] [--ablation]
Resumable; writes results/raw_runs.csv, fronts.csv, hv_history.csv, validation.csv."""
import argparse
import csv
import json
import os
import platform
import sys
import time
from multiprocessing import Pool

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from model import (Problem, agg_fnd_true_distance, exact_agg_degree, exact_c, lifetime_agg,   # noqa: E402
                   lifetime_from_c, tree_fnd_true_distance)
from moea import nsga2_sa, pymoo_run, de_cr, eps_de, hv2d                           # noqa: E402

# A1_CONFIG / A1_RES select the /boucle configuration and output folder (defaults: main campaign)
CFG = json.load(open(os.path.join(ROOT, "experiments", os.environ.get("A1_CONFIG", "config.json"))))
S = CFG["seeds"]
RES = os.path.join(ROOT, os.environ.get("A1_RES", "results"))
CROSS = bool(CFG.get("cross_regime", False))   # store each front point's lifetime under the other regime
FIELDS = ["scenario", "N", "L", "Rc", "algorithm", "run", "target_seed", "init_seed", "algo_seed",
          "budget", "n_evals", "HV", "max_CR", "max_Tn", "n_front", "T_max_rounds", "time_s"]


def solve(alg, prob, pop0, r):
    B = CFG["budget"]
    rng = np.random.default_rng(S["algo_seed"][r])
    if alg == "NSGA-II+SA":
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"])
    if alg.startswith("NSGA-II+SA p_sa="):
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["ablation"][alg])
    if alg == "NSGA-II":
        return nsga2_sa(prob, B, pop0, rng, p_sa=0.0)
    if alg.startswith("NSGA-II ec="):     # tuned variation parameters, e.g. "NSGA-II ec=5 em=50"
        ec, em = (float(t.split("=")[1]) for t in alg.split()[1:])
        return nsga2_sa(prob, B, pop0, rng, p_sa=0.0, eta_c=ec, eta_m=em)
    if alg == "NSGA-II+RND":
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"], mode=1)
    if alg == "NSGA-II+REPAIR":
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"], mode=2)
    if alg in ("NSGA-II+KY-RH", "NSGA-II+DD-MCSA"):
        from repairs_pub import dd_mcsa, ky_rh
        f = ky_rh if alg.endswith("KY-RH") else dd_mcsa
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"], repair=lambda xy, g: f(xy, prob.sink, prob.rc, g))
    if alg == "EPS-DE":
        return eps_de(prob, B, pop0, rng)
    if alg == "EPS-DE+REPAIR":
        return eps_de(prob, B, pop0, rng, p_repair=CFG["p_sa"])
    if alg.startswith("NSGA-II+SA d="):
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"], step_max=float(alg.split("=")[1]))
    if alg.startswith("NSGA-II+SA m="):
        return nsga2_sa(prob, B, pop0, rng, p_sa=CFG["p_sa"], margin=float(alg.split("=")[1]))
    if alg.endswith("+SA"):
        return pymoo_run(alg[:-3], prob, B, pop0, S["algo_seed"][r], p_sa=CFG["p_sa"])
    if alg == "DE-CR":
        return de_cr(prob, B, pop0, rng)
    return pymoo_run(alg, prob, B, pop0, S["algo_seed"][r])


def task(args):
    key, alg, r = args
    s = CFG["scenarios"][key]
    prob = Problem(s["N"], s["T"], s["L"], CFG["rs"], s["Rc"], target_seed=S["target_seed"][r],
                   targets=s.get("targets", "uniform"), sink=s.get("sink", "centre"),
                   regime=s.get("regime", "raw"))
    pop0 = np.random.default_rng(S["init_seed"][r]).uniform(0, s["L"], (CFG["NP"], prob.dim))
    t0 = time.perf_counter()
    res = solve(alg, prob, pop0, r)
    dt = time.perf_counter() - t0
    F, X = res["F"], res["X"]
    row = dict(scenario=key, N=s["N"], L=s["L"], Rc=s["Rc"], algorithm=alg, run=r,
               target_seed=S["target_seed"][r], init_seed=S["init_seed"][r], algo_seed=S["algo_seed"][r],
               budget=CFG["budget"], n_evals=prob.n_evals, HV=hv2d(1.0 - F), max_CR=F[:, 0].max(),
               max_Tn=F[:, 1].max(), n_front=len(F), T_max_rounds=prob.T_max, time_s=round(dt, 2))
    fronts = [(key, alg, r, f[0], f[1]) for f in F]
    if CROSS:
        # lifetime of the same deployments under the other collection regime (not counted as evaluations)
        other = Problem(s["N"], s["T"], s["L"], CFG["rs"], s["Rc"], target_seed=S["target_seed"][r],
                        targets=s.get("targets", "uniform"), sink=s.get("sink", "centre"),
                        regime="raw" if prob.regime == "agg" else "agg")
        Fo = other.evaluate(X)[0]
        fronts = [(key, alg, r, f[0], f[1], fo[1]) for f, fo in zip(F, Fo)]
    hist = [(key, alg, r, int(a), b) for a, b in res["hist"]]
    val = []
    pick = np.unique(np.linspace(0, len(F) - 1, min(10, len(F))).astype(int))
    order = np.argsort(F[:, 0])
    for i in order[pick]:
        xy = X[i].reshape(-1, 2)
        if F[i, 1] <= 0:
            continue
        T_tree = F[i, 1] * prob.T_max
        if prob.regime == "agg":
            val.append((key, alg, r, F[i, 0], T_tree, lifetime_agg(exact_agg_degree(xy, prob.sink, prob.rc), prob.rc),
                        agg_fnd_true_distance(xy, prob.sink, prob.rc)))
            continue
        c_star = exact_c(xy, prob.sink, prob.rc)
        val.append((key, alg, r, F[i, 0], T_tree, lifetime_from_c(c_star, prob.rc),
                    tree_fnd_true_distance(xy, prob.sink, prob.rc)))
    return row, fronts, hist, val


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=os.cpu_count())
    ap.add_argument("--ablation", action="store_true")
    ap.add_argument("--exp")
    a = ap.parse_args()
    os.makedirs(RES, exist_ok=True)
    raw = os.path.join(RES, "raw_runs.csv")
    done = set()
    if os.path.exists(raw):
        with open(raw) as fh:
            done = {(d["scenario"], d["algorithm"], int(d["run"])) for d in csv.DictReader(fh)}
    algs = list(CFG["ablation"]) if a.ablation else CFG["algorithms"]
    # slowest algorithm (MOEA/D) first for load balancing
    algs = sorted(algs, key=lambda x: x != "MOEA/D")
    keys = list(CFG["scenarios"])
    if a.exp:
        algs = sorted(CFG["algorithms_" + a.exp], key=lambda x: x != "MOEA/D")
        keys = list(CFG["scenarios_" + a.exp])
    todo = [(k, al, r) for al in algs for k in keys for r in range(CFG["runs"])
            if (k, al, r) not in done]
    import numba, scipy, pymoo
    json.dump(dict(python=sys.version, platform=platform.platform(), numpy=np.__version__,
                   numba=numba.__version__, scipy=scipy.__version__, pymoo=pymoo.__version__,
                   started=time.strftime("%Y-%m-%d %H:%M:%S")),
              open(os.path.join(RES, "environment.json"), "w"), indent=1)
    print(f"{len(todo)} runs", flush=True)
    files = {n: os.path.join(RES, n) for n in ["fronts.csv", "hv_history.csv", "validation.csv"]}
    heads = {"fronts.csv": ["scenario", "algorithm", "run", "CR", "Tn"] + (["Tn_other"] if CROSS else []),
             "hv_history.csv": ["scenario", "algorithm", "run", "fe", "HV"],
             "validation.csv": ["scenario", "algorithm", "run", "CR", "T_tree", "T_exact_flow", "FND_true_dist"]}
    new = {n: not os.path.exists(p) for n, p in files.items()}
    new_raw = not os.path.exists(raw)
    fh = {n: open(p, "a", newline="") for n, p in files.items()}
    wr = {n: csv.writer(f) for n, f in fh.items()}
    for n in files:
        if new[n]:
            wr[n].writerow(heads[n])
    with open(raw, "a", newline="") as fr, Pool(a.workers) as pool:
        w = csv.DictWriter(fr, FIELDS)
        if new_raw:
            w.writeheader()
        for k, (row, fronts, hist, val) in enumerate(pool.imap_unordered(task, todo), 1):
            w.writerow(row)
            wr["fronts.csv"].writerows(fronts)
            wr["hv_history.csv"].writerows(hist)
            wr["validation.csv"].writerows(val)
            fr.flush()
            for f in fh.values():
                f.flush()
            if k % 20 == 0:
                print(f"{k}/{len(todo)} {time.strftime('%H:%M:%S')}", flush=True)


if __name__ == "__main__":
    main()
