### Table 6. Lifetime model validation on final-front points

| scenario      |    n |   tree_over_flow_median |   tree_over_flow_min |   fnd_true_over_tree_median |   fnd_true_over_tree_min |
|:--------------|-----:|------------------------:|---------------------:|----------------------------:|-------------------------:|
| N30_Rc12_L100 | 2544 |                       1 |                0.512 |                       1     |                        1 |
| N30_Rc12_L140 | 2522 |                       1 |                0.57  |                       1.001 |                        1 |
| N30_Rc8_L100  | 2543 |                       1 |                0.63  |                       1     |                        1 |
| N30_Rc8_L140  | 2533 |                       1 |                0.524 |                       1     |                        1 |
| N50_Rc12_L100 | 2563 |                       1 |                0.552 |                       1     |                        1 |
| N50_Rc12_L140 | 2522 |                       1 |                0.512 |                       1     |                        1 |
| N50_Rc8_L100  | 2532 |                       1 |                0.548 |                       1     |                        1 |
| N50_Rc8_L140  | 2528 |                       1 |                0.387 |                       1     |                        1 |

T_tree / T_flow ≤ 1 by construction (load-balanced tree vs optimal flow routing). FND with true distances / floor(T_tree) ≥ 1 expected (all links ≤ Rc). Spearman(T_tree, FND_true) = 0.9907 over 20287 points.
