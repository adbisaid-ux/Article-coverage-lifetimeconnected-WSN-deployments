"""A1 model: bi-objective connected target coverage vs. network lifetime.

Network: N sensors (decision variables: 2N coordinates in [0, L]^2), T targets,
a sink at the field centre. Links (sensor-sensor and sensor-sink) exist when the
distance is <= Rc. Only sensors with a multi-hop path to the sink are "connected".

Objective 1 - connected coverage CR: fraction of targets within rs of a connected sensor.

Objective 2 - lifetime T under maximum-lifetime routing (flow routing, no aggregation,
one k-bit packet per connected sensor per round). Every link is <= Rc, so a
conservative per-packet cost is used for every hop:
    e_tx = k (E_elec + eps_fs Rc^2),   e_rx = k E_elec.
If c is the largest number of packets a sensor transmits per round (its own + relayed),
its energy per round is c e_tx + (c - 1) e_rx, hence
    T(c) = E0 / (c e_tx + (c - 1) e_rx),   T_max = T(1) = E0 / e_tx.
Objective used during optimisation: c_tree = largest node load of a load-balanced
minimum-hop routing tree (balanced_tree_load, O(N^2)), i.e. the lifetime of a concrete,
implementable routing. Reference: c* = smallest uniform node capacity for which all
packets reach the sink with unrestricted (flow) routing, computed by max-flow and binary
search on the final solutions only; c* <= c_tree, so T(c_tree) <= T(c*), and the gap
measures what an optimal routing could still add.
"""
import numpy as np
from numba import njit
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_flow

RADIO = dict(E0=0.5, k_bits=4000, E_elec=50e-9, eps_fs=10e-12)


def radio_costs(rc, radio=RADIO):
    e_tx = radio["k_bits"] * (radio["E_elec"] + radio["eps_fs"] * rc * rc)
    e_rx = radio["k_bits"] * radio["E_elec"]
    return e_tx, e_rx


def lifetime_from_c(c, rc, radio=RADIO):
    e_tx, e_rx = radio_costs(rc, radio)
    return radio["E0"] / (c * e_tx + (c - 1.0) * e_rx)


@njit(cache=True)
def bfs_hops(xy, sink, rc2):
    """Hop count from the sink (-1 if disconnected) and BFS order."""
    n = xy.shape[0]
    hop = np.full(n, -1)
    order = np.empty(n, np.int64)
    head = 0
    tail = 0
    for i in range(n):
        dx = xy[i, 0] - sink[0]
        dy = xy[i, 1] - sink[1]
        if dx * dx + dy * dy <= rc2:
            hop[i] = 1
            order[tail] = i
            tail += 1
    while head < tail:
        u = order[head]
        head += 1
        for v in range(n):
            if hop[v] < 0:
                dx = xy[u, 0] - xy[v, 0]
                dy = xy[u, 1] - xy[v, 1]
                if dx * dx + dy * dy <= rc2:
                    hop[v] = hop[u] + 1
                    order[tail] = v
                    tail += 1
    return hop, order[:tail]


@njit(cache=True)
def balanced_tree_load(xy, hop, order, rc2):
    """Load-balanced minimum-hop tree: nodes are processed from the farthest hop inwards
    and each picks, among its neighbours one hop closer, the one with the smallest
    current load (own packet + relayed packets). Returns the largest node load, i.e. the
    capacity c_tree of a concrete routing tree (c* <= c_tree)."""
    n = xy.shape[0]
    m = order.shape[0]
    load = np.zeros(n)
    for k in range(m):
        load[order[k]] = 1.0
    for k in range(m - 1, -1, -1):
        v = order[k]
        if hop[v] == 1:
            continue
        best = 1e300
        bu = -1
        for u in range(n):
            if hop[u] == hop[v] - 1:
                dx = xy[u, 0] - xy[v, 0]
                dy = xy[u, 1] - xy[v, 1]
                if dx * dx + dy * dy <= rc2 and load[u] < best:
                    best = load[u]
                    bu = u
        load[bu] += load[v]
    return load.max()


@njit(cache=True)
def evaluate_one(xy, targets, sink, rs2, rc2, agg=False):
    """Returns CR, c_tree (or d_max in the aggregation regime; 0 if nothing connected),
    number of connected sensors."""
    hop, order = bfs_hops(xy, sink, rc2)
    m = order.shape[0]
    if m == 0:
        return 0.0, 0.0, 0
    cov = 0
    for t in range(targets.shape[0]):
        for i in range(xy.shape[0]):
            if hop[i] > 0:
                dx = xy[i, 0] - targets[t, 0]
                dy = xy[i, 1] - targets[t, 1]
                if dx * dx + dy * dy <= rs2:
                    cov += 1
                    break
    if agg:
        return cov / targets.shape[0], agg_tree_degree(xy, hop, order, rc2), m
    c_tree = balanced_tree_load(xy, hop, order, rc2)
    return cov / targets.shape[0], c_tree, m


@njit(cache=True)
def evaluate_batch(pop, targets, sink, rs, rc, agg=False):
    p = pop.shape[0]
    n = pop.shape[1] // 2
    cr = np.empty(p)
    clb = np.empty(p)
    m = np.empty(p, np.int64)
    for i in range(p):
        a, b, c = evaluate_one(pop[i].reshape((n, 2)), targets, sink, rs * rs, rc * rc, agg)
        cr[i] = a
        clb[i] = b
        m[i] = c
    return cr, clb, m


class Problem:
    def __init__(self, n_sensors, n_targets, area, rs, rc, target_seed, radio=RADIO,
                 targets="uniform", sink="centre", regime="raw"):
        self.N, self.T, self.L, self.rs, self.rc = n_sensors, n_targets, float(area), float(rs), float(rc)
        self.dim = 2 * n_sensors
        self.lower, self.upper = 0.0, float(area)
        rng = np.random.default_rng(target_seed)
        if targets == "uniform":
            self.targets = rng.uniform(0, area, (n_targets, 2))
        elif targets == "clustered":      # 5 Gaussian clusters, sd = area/10, centres in the inner 80 %
            cen = rng.uniform(0.1 * area, 0.9 * area, (5, 2))
            lab = rng.integers(5, size=n_targets)
            self.targets = np.clip(cen[lab] + rng.normal(0, area / 10, (n_targets, 2)), 0, area)
        else:
            raise ValueError(targets)
        self.sink = {"centre": np.array([area / 2.0, area / 2.0]), "corner": np.array([0.0, 0.0])}[sink]
        self.radio = radio
        self.regime = regime          # "raw": no aggregation (c_tree); "agg": aggregation (d_max)
        self.T_max = lifetime_from_c(1.0, rc, radio)
        self.n_evals = 0

    def evaluate(self, pop):
        """Objectives to MAXIMISE: (CR, T/T_max). T = 0 when nothing is connected."""
        pop = np.ascontiguousarray(np.atleast_2d(pop), dtype=np.float64)
        agg = self.regime == "agg"
        cr, clb, m = evaluate_batch(pop, self.targets, self.sink, self.rs, self.rc, agg)
        if agg:
            t = np.where(m > 0, lifetime_agg(clb, self.rc, self.radio), 0.0)
        else:
            t = np.where(m > 0, lifetime_from_c(np.maximum(clb, 1.0), self.rc, self.radio), 0.0)
        self.n_evals += pop.shape[0]
        return np.column_stack([cr, t / self.T_max]), clb, m


# ---------------------------------------------------------------- exact value (validation)
def exact_c(xy, sink, rc, scale=1000, tol=1e-3):
    """Smallest uniform node capacity c* (packets/round) routing every connected node's
    packet to the sink, by binary search on c with integer max-flow (node splitting)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    hop, order = bfs_hops(xy, sink, rc * rc)
    conn = np.flatnonzero(hop > 0)
    m = conn.size
    if m == 0:
        return np.nan
    idx = {v: k for k, v in enumerate(conn)}
    # nodes: S=0, in_k = 1+k, out_k = 1+m+k, sink = 1+2m
    S, SINK = 0, 1 + 2 * m
    d2 = ((xy[conn, None, :] - xy[None, conn, :]) ** 2).sum(-1)
    adj = (d2 <= rc * rc) & ~np.eye(m, dtype=bool)
    near_sink = (((xy[conn] - sink) ** 2).sum(-1) <= rc * rc)

    def feasible(c):
        rows, cols, caps = [], [], []
        cap_node = int(np.floor(c * scale + 1e-9))
        for k in range(m):
            rows += [S, 1 + k]; cols += [1 + k, 1 + m + k]; caps += [scale, cap_node]
            if near_sink[k]:
                rows.append(1 + m + k); cols.append(SINK); caps.append(m * scale)
            for j in np.flatnonzero(adj[k]):
                rows.append(1 + m + k); cols.append(1 + j); caps.append(m * scale)
        g = csr_matrix((np.array(caps, np.int32), (rows, cols)), shape=(2 + 2 * m, 2 + 2 * m))
        return maximum_flow(g, S, SINK).flow_value >= m * scale

    lo, hi = 1.0, float(m)
    if feasible(lo):
        return 1.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if feasible(mid):
            hi = mid
        else:
            lo = mid
    return hi


# ---------------------------------------------------------------- sensor-agent operator
@njit(cache=True)
def agent_step(xy_in, targets, sink, rs, rc, lo, hi, step_max, rnd, mode=0, margin=0.95):
    """SA-op: every sensor is an agent applying a fixed rule; the tree and the uncovered targets are
    computed centrally for the whole deployment. Agents are processed from the
    sink outwards (BFS order), so a parent always moves before its children.
      - redundant connected agent (no target covered only by itself): moves towards the
        nearest uncovered target, at most step_max, stopping once the target is sensed;
      - every connected agent is then projected into the disc of radius 0.95 Rc around
        its parent's final position (the sink for hop-1 agents): connectivity of the BFS
        tree is preserved by construction;
      - disconnected agents finally attach to the nearest connected agent, at 0.9 Rc in
        the direction of their current position (random direction if co-located).
    rnd: array of uniform numbers (>= 2N) for random directions.
    mode (control variants for the ablation, /boucle 17/10/2026):
      0 = SA-op as above; 1 = redundant agents move by step_max in a random direction
      (same projection and attachment); 2 = no move, projection and attachment only (repair)."""
    # attachment radius: rounded so that margin = 0.95 gives exactly 0.9, as in the main campaign
    att = round(margin - 0.05, 10)
    xy = xy_in.copy()
    n = xy.shape[0]
    rs2 = rs * rs
    rc2 = rc * rc
    hop, order = bfs_hops(xy, sink, rc2)
    # parent = nearest neighbour one hop closer (-1 = sink)
    parent = np.full(n, -1)
    for v in order:
        if hop[v] > 1:
            best = 1e300
            for u in range(n):
                if hop[u] == hop[v] - 1:
                    dx = xy[u, 0] - xy[v, 0]
                    dy = xy[u, 1] - xy[v, 1]
                    d = dx * dx + dy * dy
                    if d <= rc2 and d < best:
                        best = d
                        parent[v] = u
    nt = targets.shape[0]
    count = np.zeros(nt, np.int64)
    for t in range(nt):
        for i in range(n):
            if hop[i] > 0:
                dx = xy[i, 0] - targets[t, 0]
                dy = xy[i, 1] - targets[t, 1]
                if dx * dx + dy * dy <= rs2:
                    count[t] += 1
    for v in order:
        # does v cover some target alone?
        unique = False
        for t in range(nt):
            if count[t] == 1:
                dx = xy[v, 0] - targets[t, 0]
                dy = xy[v, 1] - targets[t, 1]
                if dx * dx + dy * dy <= rs2:
                    unique = True
                    break
        old = xy[v].copy()
        if mode == 1 and not unique:
            ang = 2 * np.pi * rnd[(n + v) % rnd.shape[0]]
            xy[v, 0] += step_max * np.cos(ang)
            xy[v, 1] += step_max * np.sin(ang)
        elif mode == 0 and not unique:
            best = 1e300
            bt = -1
            for t in range(nt):
                if count[t] == 0:
                    dx = targets[t, 0] - xy[v, 0]
                    dy = targets[t, 1] - xy[v, 1]
                    d = dx * dx + dy * dy
                    if d < best:
                        best = d
                        bt = t
            if bt >= 0:
                d = np.sqrt(best)
                s = min(step_max, max(d - 0.9 * rs, 0.0))
                if d > 0:
                    xy[v, 0] += s * (targets[bt, 0] - xy[v, 0]) / d
                    xy[v, 1] += s * (targets[bt, 1] - xy[v, 1]) / d
        # projection onto the parent's disc
        if parent[v] >= 0:
            px, py = xy[parent[v], 0], xy[parent[v], 1]
        else:
            px, py = sink[0], sink[1]
        dx = xy[v, 0] - px
        dy = xy[v, 1] - py
        d = np.sqrt(dx * dx + dy * dy)
        r = margin * rc
        if d > r:
            xy[v, 0] = px + dx * r / d
            xy[v, 1] = py + dy * r / d
        xy[v, 0] = min(max(xy[v, 0], lo), hi)
        xy[v, 1] = min(max(xy[v, 1], lo), hi)
        # bounds may break the parent link only if the parent disc leaves the field; re-check
        dx = xy[v, 0] - px
        dy = xy[v, 1] - py
        if dx * dx + dy * dy > rc2:
            xy[v, 0] = px
            xy[v, 1] = py
        # update coverage counts
        for t in range(nt):
            dxo = old[0] - targets[t, 0]
            dyo = old[1] - targets[t, 1]
            dxn = xy[v, 0] - targets[t, 0]
            dyn = xy[v, 1] - targets[t, 1]
            a = dxo * dxo + dyo * dyo <= rs2
            b = dxn * dxn + dyn * dyn <= rs2
            if a and not b:
                count[t] -= 1
            elif b and not a:
                count[t] += 1
    # attach disconnected agents
    connected = hop > 0
    k = 0
    for v in range(n):
        if not connected[v]:
            best = 1e300
            bu = -1
            for u in range(n):
                if connected[u]:
                    dx = xy[u, 0] - xy[v, 0]
                    dy = xy[u, 1] - xy[v, 1]
                    d = dx * dx + dy * dy
                    if d < best:
                        best = d
                        bu = u
            if bu < 0:    # nothing connected yet: attach to the sink
                ax, ay = sink[0], sink[1]
                dx = xy[v, 0] - ax
                dy = xy[v, 1] - ay
            else:
                ax, ay = xy[bu, 0], xy[bu, 1]
                dx = xy[v, 0] - ax
                dy = xy[v, 1] - ay
            d = np.sqrt(dx * dx + dy * dy)
            if d < 1e-9:
                ang = 2 * np.pi * rnd[k % rnd.shape[0]]
                dx, dy, d = np.cos(ang), np.sin(ang), 1.0
            k += 1
            nx = min(max(ax + att * rc * dx / d, lo), hi)
            ny = min(max(ay + att * rc * dy / d, lo), hi)
            ddx = nx - ax
            ddy = ny - ay
            if ddx * ddx + ddy * ddy > rc2:
                nx, ny = ax, ay
            xy[v, 0] = nx
            xy[v, 1] = ny
            connected[v] = True
    return xy


def tree_fnd_true_distance(xy, sink, rc, radio=RADIO):
    """FND (rounds) of the load-balanced tree with true link distances (validation)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    n = xy.shape[0]
    hop, order = bfs_hops(xy, sink, rc * rc)
    if order.size == 0:
        return np.nan
    load = np.zeros(n)
    load[order] = 1.0
    parent = np.full(n, -1)
    for v in order[::-1]:
        if hop[v] == 1:
            continue
        cand = np.flatnonzero((hop == hop[v] - 1) & (((xy - xy[v]) ** 2).sum(1) <= rc * rc))
        u = cand[np.argmin(load[cand])]       # same rule and tie-breaking as balanced_tree_load
        load[u] += load[v]
        parent[v] = u
    k, Ee, efs = radio["k_bits"], radio["E_elec"], radio["eps_fs"]
    e = np.zeros(n)
    for v in order:
        d = np.hypot(*(xy[v] - (sink if parent[v] < 0 else xy[parent[v]])))
        e[v] = load[v] * k * (Ee + efs * d * d) + (load[v] - 1) * k * Ee
    return float(np.floor(radio["E0"] / e[order].max()))


# ---------------------------------------------------------------- aggregation regime
@njit(cache=True)
def agg_tree_assign(xy, hop, order, rc2):
    """Aggregation regime (cluster-head style data fusion): every relay merges the packets of
    its children with its own reading and transmits ONE packet per round. Its energy per round
    is e_tx + d e_rx, where d is its number of children, so lifetime is set by the largest
    in-degree d_max of the routing tree. Among minimum-hop trees, d_max is minimised exactly:
    each sensor at hop h >= 2 must choose a parent at hop h-1 within Rc, a semi-matching
    problem solved optimally by the incremental alternating-path algorithm of Harvey et al.
    (each new sensor is assigned to the least-loaded parent reachable by an alternating path).
    The sink has no energy limit, so hop-1 sensors are not counted as children of anyone.
    Returns the number of children of every sensor and the parent of every sensor at hop >= 2
    (-1 otherwise: hop-1 sensors send to the sink)."""
    n = xy.shape[0]
    m = order.shape[0]
    load = np.zeros(n, np.int64)
    assign = np.full(n, -1)
    via = np.full(n, -1)
    seen = np.zeros(n, np.int64)
    queue = np.empty(n, np.int64)
    stamp = 0
    for k in range(m):
        v = order[k]
        if hop[v] < 2:
            continue
        stamp += 1
        head = 0
        tail = 0
        best = -1
        for u in range(n):
            if hop[u] == hop[v] - 1:
                dx = xy[u, 0] - xy[v, 0]
                dy = xy[u, 1] - xy[v, 1]
                if dx * dx + dy * dy <= rc2:
                    seen[u] = stamp
                    via[u] = v
                    queue[tail] = u
                    tail += 1
        while head < tail:
            u = queue[head]
            head += 1
            if best < 0 or load[u] < load[best]:
                best = u
            for w in range(n):
                if assign[w] == u:
                    for u2 in range(n):
                        if hop[u2] == hop[u] and seen[u2] != stamp:
                            dx = xy[u2, 0] - xy[w, 0]
                            dy = xy[u2, 1] - xy[w, 1]
                            if dx * dx + dy * dy <= rc2:
                                seen[u2] = stamp
                                via[u2] = w
                                queue[tail] = u2
                                tail += 1
        load[best] += 1
        u = best
        while True:
            t = via[u]
            old = assign[t]
            assign[t] = u
            if t == v:
                break
            u = old
    return load, assign


@njit(cache=True)
def agg_tree_degree(xy, hop, order, rc2):
    """Largest in-degree of the optimal aggregation tree (see agg_tree_assign)."""
    if order.shape[0] == 0:
        return 0.0
    load, assign = agg_tree_assign(xy, hop, order, rc2)
    return float(load.max())


def agg_fnd_true_distance(xy, sink, rc, radio=RADIO):
    """First-node-death lifetime (rounds) of the aggregation tree with every link charged at
    its true length instead of Rc (>= the model lifetime)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    hop, order = bfs_hops(xy, sink, rc * rc)
    if order.size == 0:
        return np.nan
    load, assign = agg_tree_assign(xy, hop, order, rc * rc)
    k, Ee, efs = radio["k_bits"], radio["E_elec"], radio["eps_fs"]
    worst = 0.0
    for v in order:
        p = sink if assign[v] < 0 else xy[assign[v]]
        d2 = ((xy[v] - p) ** 2).sum()
        worst = max(worst, k * (Ee + efs * d2) + load[v] * k * Ee)
    return radio["E0"] / worst


def lifetime_agg(d, rc, radio=RADIO):
    e_tx, e_rx = radio_costs(rc, radio)
    return radio["E0"] / (e_tx + d * e_rx)


def exact_agg_degree(xy, sink, rc):
    """Validation of agg_tree_degree: smallest D such that every sensor at hop >= 2 can be
    given a parent one hop closer with at most D children per parent (max-flow, D = 0, 1, ...)."""
    xy = np.asarray(xy, float).reshape(-1, 2)
    hop, order = bfs_hops(xy, sink, rc * rc)
    tasks = np.flatnonzero(hop >= 2)
    if tasks.size == 0:
        return 0
    mach = np.flatnonzero(hop >= 1)
    mi = {u: k for k, u in enumerate(mach)}
    nt, nm = tasks.size, mach.size
    S, T = 0, 1 + nt + nm
    for D in range(1, nt + 1):
        rows, cols, caps = [], [], []
        for a, v in enumerate(tasks):
            rows.append(S); cols.append(1 + a); caps.append(1)
            for u in mach[hop[mach] == hop[v] - 1]:
                if ((xy[u] - xy[v]) ** 2).sum() <= rc * rc:
                    rows.append(1 + a); cols.append(1 + nt + mi[u]); caps.append(1)
        for b in range(nm):
            rows.append(1 + nt + b); cols.append(T); caps.append(D)
        g = csr_matrix((np.array(caps, np.int32), (rows, cols)), shape=(T + 1, T + 1))
        if maximum_flow(g, S, T).flow_value == nt:
            return D
    return nt
