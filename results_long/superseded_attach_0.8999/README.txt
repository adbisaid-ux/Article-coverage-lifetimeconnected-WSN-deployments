Runs made between 17/10 15:18 and 18/10 00:10 with agent_step attaching disconnected sensors at (0.95-0.05)*Rc = 0.8999999999999999*Rc
instead of 0.9*Rc as in the main campaign (floating-point side effect of the margin parameter). Fixed on 18/10/2026
(att = round(margin - 0.05, 10)); the affected runs (SA-op and REPAIR with margin 0.95) were re-run with the fixed code.
These copies are kept for the record only.

Comparison after re-running (18/10/2026): the REPAIR-based runs (EPS-DE+REPAIR in R3, NSGA-II+REPAIR in L) came out
bit-for-bit identical, so the rounding did not change their trajectories; only the SA-op runs (R2 d=4, d=16 and L
NSGA-II+SA) changed. The one conclusion that moved is the L comparison NSGA-II+SA vs NSGA-II in N50_Rc12_L100
(adjusted p 0.31 before, 0.039 after), reported as borderline in the manuscript.
