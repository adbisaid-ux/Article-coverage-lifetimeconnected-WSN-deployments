### Front metrics averaged over the 8 scenarios x 30 runs

| algorithm           |   n_distinct |   IGDp |   HV_relay |   share_c1 |
|:--------------------|-------------:|-------:|-----------:|-----------:|
| NSGA-II+SA          |       9.925  | 0.0165 |     0.0939 |     0.3166 |
| NSGA-II             |       7.6667 | 0.0845 |     0.08   |     0.3652 |
| NSGA-II(pymoo)      |       6.9875 | 0.114  |     0.0726 |     0.3899 |
| NSGA-III            |       8      | 0.0913 |     0.0789 |     0.3626 |
| MOEA/D              |       5.7375 | 0.1829 |     0.0544 |     0.4893 |
| SMS-EMOA            |       6.5667 | 0.1653 |     0.0595 |     0.4372 |
| DE-CR               |       3.8792 | 0.2477 |     0.02   |     0      |
| NSGA-II+SA p_sa=0.1 |       9.7458 | 0.0183 |     0.0932 |     0.3175 |
| NSGA-II+SA p_sa=0.5 |       9.8208 | 0.0152 |     0.0935 |     0.3151 |
| NSGA-II+RND         |       9.5667 | 0.0227 |     0.0932 |     0.3158 |
| NSGA-II+REPAIR      |       9.8875 | 0.0188 |     0.0923 |     0.3184 |
| NSGA-III+SA         |      10.6292 | 0.0337 |     0.0877 |     0.3326 |
| SMS-EMOA+SA         |      10.925  | 0.037  |     0.0803 |     0.3512 |

n_distinct: distinct objective vectors on the final front; IGD+: w.r.t. the non-dominated union of all algorithms on the same instance (lower is better); HV_relay: HV of the front points with load c >= 2; share_c1: share of HV contributed by the no-relay point.
