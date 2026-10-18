"""Published connectivity repairs, re-implemented as baselines for SA-op (review round 4, item 7).

KY-RH: repair heuristic of Konstantinidis and Yang (Computer Communications, 2011, doi 10.1016/j.comcom.2010.08.011),
       adapted to K = 1, one fixed communication range Rc (no power assignment) and NSGA-II (no subproblem weight):
       - sensor j is feasible if it reaches the sink directly or has a neighbour strictly closer to the sink
         (positive advance); with K = 1 all sensors feasible implies every sensor is connected to the sink;
       - while an infeasible sensor exists, one of them (uniformly at random) is repaired:
         with probability 1/2 (the paper uses the subproblem weight k_i >= 0.5) it is moved to a uniformly random
         point of the sparsest of s = 4 sectors of the disc of radius Rc around the sink;
         otherwise (k_i < 0.5) it is moved along the line towards u, to distance Rc from u (Eq. 14 of the paper),
         where u is its nearest sensor v if v is feasible and closer to the sink than j, and the sink otherwise;
       - the loop stops when all sensors are feasible or after 20 N repairs.
DD-MCSA: connectivity phase of Deghbouch and Debbat (Informatica 46(8), 2022, doi 10.31449/inf.v46i8.3370):
       the set of connected nodes starts from the sink; each disconnected sensor moves towards a connected node
       chosen at random, x <- x + fl (m - x), with fl growing linearly from 0.1 to 1 over MaxIter = 10 iterations;
       connected sensors do not move; the loop ends when the network is connected (at the latest when fl = 1).
Both act on one deployment (N x 2 array) and return a new array; the caller applies them to offspring at the same
rate p_sa as SA-op."""
import numpy as np

from model import bfs_hops


def _feasible(xy, sink, rc2):
    """Positive-advance feasibility of every sensor (K = 1)."""
    dH = ((xy - sink) ** 2).sum(1)
    d = ((xy[:, None, :] - xy[None, :, :]) ** 2).sum(-1)
    closer = dH[None, :] < dH[:, None]                     # closer[j, i]: i is closer to the sink than j
    return (dH <= rc2) | ((d <= rc2) & closer).any(1), dH, d


def ky_rh(xy_in, sink, rc, rng, sectors=4, max_rep_per_sensor=20):
    xy = xy_in.copy()
    n = len(xy)
    rc2 = rc * rc
    for _ in range(max_rep_per_sensor * n):
        ok, dH, d = _feasible(xy, sink, rc2)
        bad = np.flatnonzero(~ok)
        if len(bad) == 0:
            break
        j = bad[rng.integers(len(bad))]
        if rng.random() < 0.5:
            # sparsest sector of the disc of radius Rc around the sink (sensors counted inside the disc)
            ang = np.arctan2(xy[:, 1] - sink[1], xy[:, 0] - sink[0]) % (2 * np.pi)
            inside = dH <= rc2
            cnt = np.bincount((ang[inside] // (2 * np.pi / sectors)).astype(int).clip(0, sectors - 1),
                              minlength=sectors)
            s = rng.choice(np.flatnonzero(cnt == cnt.min()))
            a = (s + rng.random()) * 2 * np.pi / sectors
            r = rc * np.sqrt(rng.random())                  # uniform in the sector
            xy[j] = sink + r * np.array([np.cos(a), np.sin(a)])
        else:
            dj = d[j].copy()
            dj[j] = np.inf
            v = int(np.argmin(dj))
            u = xy[v] if (ok[v] and dH[v] < dH[j]) else sink
            duj = np.sqrt(((u - xy[j]) ** 2).sum())
            if duj > 0:
                # Eq. 14: move along the line towards u, to distance Rc from u (1e-9 margin against rounding)
                xy[j] = xy[j] + (duj - rc * (1 - 1e-9)) * (u - xy[j]) / duj
    return xy


def dd_mcsa(xy_in, sink, rc, rng, max_iter=10):
    xy = xy_in.copy()
    rc2 = rc * rc
    for it in range(1, max_iter + 1):
        hop, _ = bfs_hops(xy, sink, rc2)
        dis = np.flatnonzero(hop < 0)
        if len(dis) == 0:
            break
        con = np.flatnonzero(hop > 0)
        fl = 0.1 + (1.0 - 0.1) * it / max_iter               # 0.19 ... 1.0; Eq. 6 of the paper
        for j in dis:
            k = rng.integers(len(con) + 1)                    # the sink is in the connected list
            m = sink if k == len(con) else xy[con[k]]
            xy[j] = xy[j] + fl * (m - xy[j])
    return xy
