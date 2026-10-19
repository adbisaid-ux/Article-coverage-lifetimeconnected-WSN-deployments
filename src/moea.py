"""Multi-objective optimisers for A1 (all maximise (CR, T/T_max) via minimising G = 1 - F).

- nsga2_sa : our NSGA-II (Deb et al. 2002) with SBX + polynomial mutation, plus the
             sensor-agent operator SA-op applied to each offspring with probability p_sa
             (p_sa = 0 gives plain NSGA-II, same code: ablation baseline).
- pymoo_run: NSGA-II, NSGA-III, MOEA/D, SMS-EMOA from pymoo with default settings.
- de_cr    : single-objective DE (F = 0.3, CR = 0.9, tuned in Campagne 0) on CR only,
             lifetime as tie-breaker; a "coverage-only" reference point.
Every evaluation (initial population, offspring, SA-op output) is charged to the budget.
"""
import numpy as np

from model import agent_step

REF = np.array([1.0, 1.0])     # HV reference point in G = 1 - F space


def hv2d(G, ref=REF):
    """Hypervolume of a 2-D minimisation point set w.r.t. ref."""
    G = G[(G[:, 0] < ref[0]) & (G[:, 1] < ref[1])]
    if len(G) == 0:
        return 0.0
    G = G[np.argsort(G[:, 0], kind="stable")]
    hv, best2 = 0.0, ref[1]
    for g1, g2 in G:
        if g2 < best2:
            hv += (ref[0] - g1) * (best2 - g2)
            best2 = g2
    return hv


def nondominated(G):
    """Boolean mask of non-dominated rows (minimisation)."""
    n = len(G)
    dom = (G[:, None, :] <= G[None, :, :]).all(-1) & (G[:, None, :] < G[None, :, :]).any(-1)
    return ~dom.any(0)


def nd_sort(G):
    n = len(G)
    dom = (G[:, None, :] <= G[None, :, :]).all(-1) & (G[:, None, :] < G[None, :, :]).any(-1)
    n_dom = dom.sum(0)
    rank = np.full(n, -1)
    r = 0
    front = np.flatnonzero(n_dom == 0)
    while front.size:
        rank[front] = r
        n_dom = n_dom - dom[front].sum(0)
        n_dom[rank >= 0] = -1
        front = np.flatnonzero(n_dom == 0)
        r += 1
    return rank


def crowding(G, rank):
    cd = np.zeros(len(G))
    for r in np.unique(rank):
        idx = np.flatnonzero(rank == r)
        if idx.size <= 2:
            cd[idx] = np.inf
            continue
        for k in range(G.shape[1]):
            o = idx[np.argsort(G[idx, k], kind="stable")]
            span = G[o[-1], k] - G[o[0], k]
            cd[o[0]] = cd[o[-1]] = np.inf
            if span > 0:
                cd[o[1:-1]] += (G[o[2:], k] - G[o[:-2], k]) / span
    return cd


def survive(X, G, n):
    rank = nd_sort(G)
    cd = crowding(G, rank)
    order = np.lexsort((-cd, rank))[:n]
    return X[order], G[order], rank[order], cd[order]


def sbx_pm(P1, P2, lo, hi, rng, eta_c=15.0, pc=0.9, eta_m=20.0, pm=None):
    """Vectorised SBX (bounded) + polynomial mutation, two children per pair."""
    n, D = P1.shape
    pm = 1.0 / D if pm is None else pm
    C1, C2 = P1.copy(), P2.copy()
    do = (rng.random((n, 1)) < pc) & (rng.random((n, D)) < 0.5) & (np.abs(P1 - P2) > 1e-14)
    y1, y2 = np.minimum(P1, P2), np.maximum(P1, P2)
    delta = np.where(do, y2 - y1, 1.0)
    u = rng.random((n, D))

    def bq(beta):
        alpha = 2.0 - beta ** -(eta_c + 1.0)
        return np.where(u <= 1.0 / alpha, (u * alpha) ** (1.0 / (eta_c + 1.0)),
                        (1.0 / (2.0 - u * alpha)) ** (1.0 / (eta_c + 1.0)))

    c1 = 0.5 * (y1 + y2 - bq(1.0 + 2.0 * (y1 - lo) / delta) * delta)
    c2 = 0.5 * (y1 + y2 + bq(1.0 + 2.0 * (hi - y2) / delta) * delta)
    c1, c2 = np.clip(c1, lo, hi), np.clip(c2, lo, hi)
    swap = rng.random((n, D)) < 0.5
    c1, c2 = np.where(swap, c2, c1), np.where(swap, c1, c2)
    C1, C2 = np.where(do, c1, C1), np.where(do, c2, C2)
    C = np.vstack([C1, C2])
    m = rng.random(C.shape) < pm
    if m.any():
        y = C[m]
        d1, d2 = (y - lo) / (hi - lo), (hi - y) / (hi - lo)
        u = rng.random(y.size)
        mp = 1.0 / (eta_m + 1.0)
        dq = np.where(u < 0.5, (2 * u + (1 - 2 * u) * (1 - d1) ** (eta_m + 1)) ** mp - 1,
                      1 - (2 * (1 - u) + 2 * (u - 0.5) * (1 - d2) ** (eta_m + 1)) ** mp)
        C[m] = np.clip(y + dq * (hi - lo), lo, hi)
    return C


def nsga2_sa(problem, budget, pop0, rng, p_sa=0.2, step_max=None, log_every=10, mode=0, margin=0.95,
             eta_c=15.0, eta_m=20.0, repair=None):
    lo, hi = problem.lower, problem.upper
    step = problem.rs if step_max is None else step_max
    X = pop0.copy()
    NP, D = X.shape
    F = problem.evaluate(X)[0]
    G = 1.0 - F
    fe = NP
    X, G, rank, cd = survive(X, G, NP)
    hist = [(fe, hv2d(G[nondominated(G)]))]
    gen = 0
    while fe < budget:
        m = min(NP, budget - fe)
        # binary tournaments on (rank, crowding)
        a, b = rng.integers(NP, size=(2, NP))
        win = np.where((rank[a] < rank[b]) | ((rank[a] == rank[b]) & (cd[a] > cd[b])), a, b)
        C = sbx_pm(X[win[: NP // 2]], X[win[NP // 2:]], lo, hi, rng, eta_c=eta_c, eta_m=eta_m)[:m]
        if p_sa > 0:
            sel = np.flatnonzero(rng.random(m) < p_sa)
            for i in sel:
                if repair is not None:      # published repair baseline (src/repairs_pub.py)
                    C[i] = repair(C[i].reshape(-1, 2), rng).ravel()
                    continue
                C[i] = agent_step(C[i].reshape(-1, 2), problem.targets, problem.sink, problem.rs,
                                  problem.rc, lo, hi, step, rng.random(2 * problem.N), mode, margin).ravel()
        FC = problem.evaluate(C)[0]
        fe += m
        X, G, rank, cd = survive(np.vstack([X, C]), np.vstack([G, 1.0 - FC]), NP)
        gen += 1
        if gen % log_every == 0 or fe >= budget:
            hist.append((fe, hv2d(G[nondominated(G)])))
    nd = nondominated(G)
    return dict(X=X[nd], F=1.0 - G[nd], hist=np.array(hist))


AGE_FALLBACKS = {"calls": 0, "fallbacks": 0}


def pymoo_run(name, problem, budget, pop0, seed, p_sa=0.0):
    from pymoo.core.problem import Problem as PP
    from pymoo.core.callback import Callback
    from pymoo.optimize import minimize
    from pymoo.util.ref_dirs import get_reference_directions
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.algorithms.moo.nsga3 import NSGA3
    from pymoo.algorithms.moo.moead import MOEAD
    from pymoo.algorithms.moo.sms import SMSEMOA
    from pymoo.algorithms.moo.age2 import AGEMOEA2, AGEMOEA2Survival

    # pymoo 0.6.2 bug: AGE-MOEA-II's geodesic distance divides by zero when a front point coincides
    # with the ideal point (frequent here, objectives are discrete). Fallback to Euclidean distances
    # for that call only; AGE_FALLBACKS counts the calls and the fallbacks in this process
    # (experiments/count_age_fallbacks.py).
    if not getattr(AGEMOEA2Survival, "_a1_patched", False):
        _orig = AGEMOEA2Survival.pairwise_distances

        def _safe(front, p):
            AGE_FALLBACKS["calls"] += 1
            try:
                return _orig(front, p)
            except ZeroDivisionError:
                AGE_FALLBACKS["fallbacks"] += 1
                D = np.sqrt(((front[:, None, :] - front[None, :, :]) ** 2).sum(-1))
                return np.triu(D, 1) + np.triu(D, 1).T
        AGEMOEA2Survival.pairwise_distances = staticmethod(_safe)
        AGEMOEA2Survival._a1_patched = True

    class W(PP):
        def __init__(self):
            super().__init__(n_var=problem.dim, n_obj=2, xl=problem.lower, xu=problem.upper)

        def _evaluate(self, x, out, *a, **k):
            out["F"] = 1.0 - problem.evaluate(x)[0]

    class CB(Callback):
        def __init__(self):
            super().__init__()
            self.h = []
            self.g = 0

        def notify(self, algo):
            self.g += 1
            if self.g % 10 == 0:
                Gp = algo.pop.get("F")
                self.h.append((problem.n_evals, hv2d(Gp[nondominated(Gp)])))

    NP = len(pop0)
    rd = get_reference_directions("das-dennis", 2, n_partitions=NP - 1)
    kw = {}
    if p_sa > 0:
        # SA-op after the default polynomial mutation, applied to each offspring with prob. p_sa;
        # offspring are evaluated once afterwards, so the budget accounting is unchanged.
        from pymoo.core.mutation import Mutation
        from pymoo.core.population import Population
        from pymoo.operators.mutation.pm import PM
        sa_rng = np.random.default_rng(seed + 7919)

        class PMSA(Mutation):
            def __init__(self):
                super().__init__()
                self.pm = PM(eta=20)

            def _do(self, prob, X, **k):
                X = self.pm.do(prob, Population.new(X=X), **k).get("X")
                for i in np.flatnonzero(sa_rng.random(len(X)) < p_sa):
                    X[i] = agent_step(X[i].reshape(-1, 2), problem.targets, problem.sink, problem.rs,
                                      problem.rc, problem.lower, problem.upper, problem.rs,
                                      sa_rng.random(2 * problem.N), 0).ravel()
                return X

        kw["mutation"] = PMSA()
    algo = {"NSGA-II(pymoo)": lambda: NSGA2(pop_size=NP, sampling=pop0, **kw),
            "NSGA-III": lambda: NSGA3(ref_dirs=rd, pop_size=NP, sampling=pop0, **kw),
            "MOEA/D": lambda: MOEAD(ref_dirs=rd, sampling=pop0, **kw),
            "SMS-EMOA": lambda: SMSEMOA(pop_size=NP, sampling=pop0, **kw),
            "AGE-MOEA-II": lambda: AGEMOEA2(pop_size=NP, sampling=pop0, **kw),
            "NSGA-II(pymoo,nodup)": lambda: NSGA2(pop_size=NP, sampling=pop0, eliminate_duplicates=False, **kw)}[name]()
    cb = CB()
    res = minimize(W(), algo, ("n_eval", budget), seed=seed, callback=cb, verbose=False)
    Gp = res.pop.get("F")
    Xp = res.pop.get("X")
    nd = nondominated(Gp)
    cb.h.append((problem.n_evals, hv2d(Gp[nd])))
    return dict(X=Xp[nd], F=1.0 - Gp[nd], hist=np.array(cb.h))


def de_cr(problem, budget, pop0, rng, F_=0.3, CR=0.9):
    """DE/rand/1/bin on f = CR + 1e-6 * T/T_max (coverage first, lifetime tie-break)."""
    lo, hi = problem.lower, problem.upper
    X = pop0.copy()
    NP, D = X.shape

    def fit(Z):
        Fz = problem.evaluate(Z)[0]
        return Fz[:, 0] + 1e-6 * Fz[:, 1], Fz

    f, FX = fit(X)
    fe = NP
    hist = [(fe, hv2d((1.0 - FX)[nondominated(1.0 - FX)]))]
    gen = 0
    while fe < budget:
        m = min(NP, budget - fe)
        idx = np.arange(m)
        r = np.array([rng.choice(np.delete(np.arange(NP), i), 3, replace=False) for i in idx])
        V = X[r[:, 0]] + F_ * (X[r[:, 1]] - X[r[:, 2]])
        mask = rng.random((m, D)) < CR
        mask[idx, rng.integers(D, size=m)] = True
        U = np.where(mask, V, X[:m])
        bad = (U < lo) | (U > hi)
        U[bad] = rng.uniform(lo, hi, bad.sum())
        fu, FU = fit(U)
        fe += m
        acc = fu >= f[:m]
        X[:m][acc], f[:m][acc], FX[:m][acc] = U[acc], fu[acc], FU[acc]
        gen += 1
        if gen % 10 == 0 or fe >= budget:
            hist.append((fe, hv2d((1.0 - FX)[nondominated(1.0 - FX)])))
    nd = nondominated(1.0 - FX)
    return dict(X=X[nd], F=FX[nd], hist=np.array(hist))


def eps_de(problem, budget, pop0, rng, levels=10, F_=0.3, CR=0.9, p_repair=0.0):
    """epsilon-constraint baseline: for c = 1..levels, DE/rand/1/bin maximises CR subject to
    c_tree <= c (feasibility first: infeasible solutions ranked by -c_tree), with budget/levels
    evaluations per level, all starting from pop0. Returns the non-dominated union of the final
    populations; every evaluation is charged."""
    lo, hi = problem.lower, problem.upper
    NP, D = pop0.shape
    per = budget // levels
    allX, allF, hist = [], [], []
    fe = 0

    def fit(Z, c):
        Fz, ct, m = problem.evaluate(Z)
        feas = (ct <= c + 1e-9) & (m > 0)
        return np.where(feas, Fz[:, 0], -ct / problem.N - 1.0), Fz

    for c in range(1, levels + 1):
        X = pop0.copy()
        f, FX = fit(X, c)
        used = NP
        while used < per:
            m = min(NP, per - used)
            idx = np.arange(m)
            r = np.array([rng.choice(np.delete(np.arange(NP), i), 3, replace=False) for i in idx])
            V = X[r[:, 0]] + F_ * (X[r[:, 1]] - X[r[:, 2]])
            mask = rng.random((m, D)) < CR
            mask[idx, rng.integers(D, size=m)] = True
            U = np.where(mask, V, X[:m])
            bad = (U < lo) | (U > hi)
            U[bad] = rng.uniform(lo, hi, bad.sum())
            for i in np.flatnonzero(rng.random(m) < p_repair):     # optional REPAIR (mode 2)
                U[i] = agent_step(U[i].reshape(-1, 2), problem.targets, problem.sink, problem.rs, problem.rc,
                                  lo, hi, problem.rs, rng.random(2 * problem.N), 2).ravel()
            fu, FU = fit(U, c)
            used += m
            acc = fu >= f[:m]
            X[:m][acc], f[:m][acc], FX[:m][acc] = U[acc], fu[acc], FU[acc]
        fe += used
        allX.append(X)
        allF.append(FX)
        G = 1.0 - np.vstack(allF)
        hist.append((fe, hv2d(G[nondominated(G)])))
    X, F = np.vstack(allX), np.vstack(allF)
    nd = nondominated(1.0 - F)
    return dict(X=X[nd], F=F[nd], hist=np.array(hist))
