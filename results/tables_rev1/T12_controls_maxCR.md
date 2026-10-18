### Largest coverage among control operators (paired Wilcoxon, Holm over 24 tests), scenarios in config order

| scenario      | a           | b              |     ma |     mb |      p |   p_holm | result   |
|:--------------|:------------|:---------------|-------:|-------:|-------:|---------:|:---------|
| N30_Rc8_L100  | NSGA-II+SA  | NSGA-II+RND    | 0.472  | 0.4547 | 0.0174 |   0.1565 | =        |
| N30_Rc8_L100  | NSGA-II+SA  | NSGA-II+REPAIR | 0.472  | 0.4817 | 0.17   |   1      | =        |
| N30_Rc8_L100  | NSGA-II+RND | NSGA-II+REPAIR | 0.4547 | 0.4817 | 0.0007 |   0.011  | -        |
| N30_Rc8_L140  | NSGA-II+SA  | NSGA-II+RND    | 0.2837 | 0.2727 | 0.013  |   0.1325 | =        |
| N30_Rc8_L140  | NSGA-II+SA  | NSGA-II+REPAIR | 0.2837 | 0.2873 | 0.474  |   1      | =        |
| N30_Rc8_L140  | NSGA-II+RND | NSGA-II+REPAIR | 0.2727 | 0.2873 | 0.001  |   0.015  | -        |
| N30_Rc12_L100 | NSGA-II+SA  | NSGA-II+RND    | 0.6847 | 0.6497 | 0.0002 |   0.0032 | +        |
| N30_Rc12_L100 | NSGA-II+SA  | NSGA-II+REPAIR | 0.6847 | 0.6763 | 0.1989 |   1      | =        |
| N30_Rc12_L100 | NSGA-II+RND | NSGA-II+REPAIR | 0.6497 | 0.6763 | 0.0065 |   0.0851 | =        |
| N30_Rc12_L140 | NSGA-II+SA  | NSGA-II+RND    | 0.415  | 0.4053 | 0.1237 |   0.9895 | =        |
| N30_Rc12_L140 | NSGA-II+SA  | NSGA-II+REPAIR | 0.415  | 0.4223 | 0.1621 |   1      | =        |
| N30_Rc12_L140 | NSGA-II+RND | NSGA-II+REPAIR | 0.4053 | 0.4223 | 0.012  |   0.1325 | =        |
| N50_Rc8_L100  | NSGA-II+SA  | NSGA-II+RND    | 0.67   | 0.6337 | 0      |   0.0004 | +        |
| N50_Rc8_L100  | NSGA-II+SA  | NSGA-II+REPAIR | 0.67   | 0.671  | 0.9904 |   1      | =        |
| N50_Rc8_L100  | NSGA-II+RND | NSGA-II+REPAIR | 0.6337 | 0.671  | 0.0003 |   0.0052 | -        |
| N50_Rc8_L140  | NSGA-II+SA  | NSGA-II+RND    | 0.4113 | 0.3907 | 0.0076 |   0.0912 | =        |
| N50_Rc8_L140  | NSGA-II+SA  | NSGA-II+REPAIR | 0.4113 | 0.4137 | 0.6966 |   1      | =        |
| N50_Rc8_L140  | NSGA-II+RND | NSGA-II+REPAIR | 0.3907 | 0.4137 | 0.0003 |   0.0049 | -        |
| N50_Rc12_L100 | NSGA-II+SA  | NSGA-II+RND    | 0.9127 | 0.858  | 0      |   0.0002 | +        |
| N50_Rc12_L100 | NSGA-II+SA  | NSGA-II+REPAIR | 0.9127 | 0.891  | 0.0039 |   0.0553 | =        |
| N50_Rc12_L100 | NSGA-II+RND | NSGA-II+REPAIR | 0.858  | 0.891  | 0.0003 |   0.0051 | -        |
| N50_Rc12_L140 | NSGA-II+SA  | NSGA-II+RND    | 0.615  | 0.5643 | 0      |   0.0002 | +        |
| N50_Rc12_L140 | NSGA-II+SA  | NSGA-II+REPAIR | 0.615  | 0.602  | 0.2382 |   1      | =        |
| N50_Rc12_L140 | NSGA-II+RND | NSGA-II+REPAIR | 0.5643 | 0.602  | 0.0001 |   0.0016 | -        |

NSGA-II+SA vs NSGA-II+RND: ==+=+=++
NSGA-II+SA vs NSGA-II+REPAIR: ========
NSGA-II+RND vs NSGA-II+REPAIR: --==----
