Reproducibility: 720 XR runs matched to earlier runs; max |HV - HV_ref| = 0.000e+00

### Aggregation regime: HV mean ± sd over 30 runs

| sc            | AGE-MOEA-II     | NSGA-II         | NSGA-II+REPAIR   | NSGA-II+RND     | NSGA-II+SA      | NSGA-III        | SMS-EMOA        |
|:--------------|:----------------|:----------------|:-----------------|:----------------|:----------------|:----------------|:----------------|
| N30_Rc8_L100  | 0.1573 ± 0.0294 | 0.1956 ± 0.0257 | 0.2659 ± 0.0238  | 0.2609 ± 0.0248 | 0.2654 ± 0.0279 | 0.166 ± 0.0381  | 0.1282 ± 0.0223 |
| N30_Rc8_L140  | 0.0671 ± 0.0249 | 0.0794 ± 0.0241 | 0.158 ± 0.0244   | 0.1552 ± 0.0218 | 0.157 ± 0.0219  | 0.0924 ± 0.0274 | 0.0588 ± 0.0231 |
| N30_Rc12_L100 | 0.296 ± 0.0384  | 0.3354 ± 0.0339 | 0.3746 ± 0.0288  | 0.3768 ± 0.0226 | 0.392 ± 0.0247  | 0.2638 ± 0.0365 | 0.2421 ± 0.0299 |
| N30_Rc12_L140 | 0.1306 ± 0.0315 | 0.167 ± 0.0254  | 0.2303 ± 0.0229  | 0.2267 ± 0.0267 | 0.2371 ± 0.0267 | 0.1486 ± 0.0294 | 0.1098 ± 0.0282 |
| N50_Rc8_L100  | 0.1958 ± 0.0364 | 0.2425 ± 0.0314 | 0.3375 ± 0.0186  | 0.3284 ± 0.023  | 0.3373 ± 0.0217 | 0.1942 ± 0.0382 | 0.1583 ± 0.026  |
| N50_Rc8_L140  | 0.0766 ± 0.0273 | 0.1037 ± 0.0215 | 0.2056 ± 0.0215  | 0.199 ± 0.0216  | 0.2046 ± 0.0236 | 0.104 ± 0.0275  | 0.0643 ± 0.023  |
| N50_Rc12_L100 | 0.3883 ± 0.0359 | 0.4293 ± 0.0288 | 0.4625 ± 0.0225  | 0.4587 ± 0.0252 | 0.4676 ± 0.0287 | 0.3043 ± 0.031  | 0.3027 ± 0.027  |
| N50_Rc12_L140 | 0.1868 ± 0.0334 | 0.2293 ± 0.0313 | 0.2983 ± 0.0242  | 0.2894 ± 0.0173 | 0.2951 ± 0.0246 | 0.1645 ± 0.0303 | 0.1378 ± 0.0222 |

### Aggregation regime: largest connected coverage, mean over 30 runs

| sc            |   AGE-MOEA-II |   NSGA-II |   NSGA-II+REPAIR |   NSGA-II+RND |   NSGA-II+SA |   NSGA-III |   SMS-EMOA |
|:--------------|--------------:|----------:|-----------------:|--------------:|-------------:|-----------:|-----------:|
| N30_Rc8_L100  |         0.244 |     0.324 |            0.5   |         0.483 |        0.497 |      0.262 |      0.186 |
| N30_Rc8_L140  |         0.099 |     0.124 |            0.304 |         0.295 |        0.299 |      0.149 |      0.082 |
| N30_Rc12_L100 |         0.484 |     0.569 |            0.679 |         0.685 |        0.71  |      0.426 |      0.387 |
| N30_Rc12_L140 |         0.204 |     0.278 |            0.437 |         0.42  |        0.446 |      0.241 |      0.161 |
| N50_Rc8_L100  |         0.324 |     0.437 |            0.701 |         0.676 |        0.699 |      0.327 |      0.257 |
| N50_Rc8_L140  |         0.118 |     0.173 |            0.445 |         0.426 |        0.445 |      0.175 |      0.093 |
| N50_Rc12_L100 |         0.703 |     0.816 |            0.915 |         0.901 |        0.931 |      0.554 |      0.565 |
| N50_Rc12_L140 |         0.325 |     0.426 |            0.631 |         0.599 |        0.627 |      0.282 |      0.229 |

HV: NSGA-II+SA vs others, Holm over 48: NSGA-II: 8+/0=/0- (gain 8.9 to 97.8 %), NSGA-II+REPAIR: 1+/7=/0- (gain -1.1 to 4.7 %), NSGA-II+RND: 2+/6=/0- (gain 1.2 to 4.6 %), NSGA-III: 8+/0=/0- (gain 48.6 to 96.7 %), SMS-EMOA: 8+/0=/0- (gain 54.5 to 217.9 %), AGE-MOEA-II: 8+/0=/0- (gain 20.4 to 166.9 %)

max_CR: NSGA-II+SA vs others, Holm over 48: NSGA-II: 8+/0=/0- (gain 14.2 to 157.2 %), NSGA-II+REPAIR: 1+/7=/0- (gain -1.5 to 4.5 %), NSGA-II+RND: 6+/2=/0- (gain 1.5 to 6.2 %), NSGA-III: 8+/0=/0- (gain 66.7 to 153.8 %), SMS-EMOA: 8+/0=/0- (gain 64.9 to 376.8 %), AGE-MOEA-II: 8+/0=/0- (gain 32.5 to 276.1 %)

### Controls in the aggregation regime (Holm over 24 tests per metric)

| scenario      | metric   | a              | b              |     ma |     mb |   p_holm | res   |    A12 |   gain% |
|:--------------|:---------|:---------------|:---------------|-------:|-------:|---------:|:------|-------:|--------:|
| N30_Rc8_L100  | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.2654 | 0.2659 |   1      | =     | 0.4922 | -0.2005 |
| N30_Rc8_L140  | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.157  | 0.158  |   1      | =     | 0.4817 | -0.6403 |
| N30_Rc12_L100 | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.392  | 0.3746 |   0.0023 | +     | 0.6867 |  4.6696 |
| N30_Rc12_L140 | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.2371 | 0.2303 |   1      | =     | 0.5722 |  2.9501 |
| N50_Rc8_L100  | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.3373 | 0.3375 |   1      | =     | 0.4922 | -0.0353 |
| N50_Rc8_L140  | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.2046 | 0.2056 |   1      | =     | 0.5072 | -0.5141 |
| N50_Rc12_L100 | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.4676 | 0.4625 |   1      | =     | 0.5922 |  1.1077 |
| N50_Rc12_L140 | HV       | NSGA-II+SA     | NSGA-II+REPAIR | 0.2951 | 0.2983 |   1      | =     | 0.4422 | -1.0734 |
| N30_Rc8_L100  | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.2654 | 0.2609 |   1      | =     | 0.5606 |  1.7301 |
| N30_Rc8_L140  | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.157  | 0.1552 |   1      | =     | 0.5333 |  1.1535 |
| N30_Rc12_L100 | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.392  | 0.3768 |   0.0011 | +     | 0.7244 |  4.0471 |
| N30_Rc12_L140 | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.2371 | 0.2267 |   0.1446 | =     | 0.5989 |  4.5847 |
| N50_Rc8_L100  | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.3373 | 0.3284 |   0.0314 | +     | 0.6106 |  2.7186 |
| N50_Rc8_L140  | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.2046 | 0.199  |   0.8787 | =     | 0.59   |  2.7903 |
| N50_Rc12_L100 | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.4676 | 0.4587 |   0.1077 | =     | 0.6056 |  1.928  |
| N50_Rc12_L140 | HV       | NSGA-II+SA     | NSGA-II+RND    | 0.2951 | 0.2894 |   1      | =     | 0.5578 |  1.9573 |
| N30_Rc8_L100  | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.2659 | 0.2609 |   1      | =     | 0.5756 |  1.9345 |
| N30_Rc8_L140  | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.158  | 0.1552 |   1      | =     | 0.5511 |  1.8054 |
| N30_Rc12_L100 | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.3746 | 0.3768 |   1      | =     | 0.5089 | -0.5947 |
| N30_Rc12_L140 | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.2303 | 0.2267 |   1      | =     | 0.5606 |  1.5877 |
| N50_Rc8_L100  | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.3375 | 0.3284 |   0.0627 | =     | 0.6256 |  2.7549 |
| N50_Rc8_L140  | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.2056 | 0.199  |   0.2817 | =     | 0.6094 |  3.3214 |
| N50_Rc12_L100 | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.4625 | 0.4587 |   1      | =     | 0.5278 |  0.8112 |
| N50_Rc12_L140 | HV       | NSGA-II+REPAIR | NSGA-II+RND    | 0.2983 | 0.2894 |   0.8618 | =     | 0.6333 |  3.0636 |
| N30_Rc8_L100  | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.497  | 0.5003 |   1      | =     | 0.4572 | -0.6662 |
| N30_Rc8_L140  | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.2993 | 0.304  |   1      | =     | 0.4611 | -1.5351 |
| N30_Rc12_L100 | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.7097 | 0.6793 |   0.0207 | +     | 0.71   |  4.4652 |
| N30_Rc12_L140 | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.446  | 0.4373 |   1      | =     | 0.565  |  1.9817 |
| N50_Rc8_L100  | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.6987 | 0.701  |   1      | =     | 0.4567 | -0.3329 |
| N50_Rc8_L140  | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.445  | 0.4447 |   1      | =     | 0.505  |  0.075  |
| N50_Rc12_L100 | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.9313 | 0.915  |   0.2131 | =     | 0.69   |  1.7851 |
| N50_Rc12_L140 | max_CR   | NSGA-II+SA     | NSGA-II+REPAIR | 0.6273 | 0.631  |   1      | =     | 0.4533 | -0.5811 |
| N30_Rc8_L100  | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.497  | 0.4827 |   0.2131 | =     | 0.6011 |  2.9696 |
| N30_Rc8_L140  | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.2993 | 0.295  |   1      | =     | 0.5411 |  1.4689 |
| N30_Rc12_L100 | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.7097 | 0.685  |   0.0459 | +     | 0.6828 |  3.601  |
| N30_Rc12_L140 | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.446  | 0.42   |   0.0312 | +     | 0.695  |  6.1905 |
| N50_Rc8_L100  | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.6987 | 0.676  |   0.021  | +     | 0.64   |  3.3531 |
| N50_Rc8_L140  | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.445  | 0.4257 |   0.035  | +     | 0.6239 |  4.5419 |
| N50_Rc12_L100 | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.9313 | 0.9013 |   0.0036 | +     | 0.7672 |  3.3284 |
| N50_Rc12_L140 | max_CR   | NSGA-II+SA     | NSGA-II+RND    | 0.6273 | 0.599  |   0.0224 | +     | 0.7028 |  4.7301 |
| N30_Rc8_L100  | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.5003 | 0.4827 |   0.0944 | =     | 0.625  |  3.6602 |
| N30_Rc8_L140  | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.304  | 0.295  |   0.2276 | =     | 0.5889 |  3.0508 |
| N30_Rc12_L100 | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.6793 | 0.685  |   1      | =     | 0.475  | -0.8273 |
| N30_Rc12_L140 | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.4373 | 0.42   |   0.2276 | =     | 0.6444 |  4.127  |
| N50_Rc8_L100  | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.701  | 0.676  |   0.0131 | +     | 0.6961 |  3.6982 |
| N50_Rc8_L140  | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.4447 | 0.4257 |   0.0065 | +     | 0.6578 |  4.4636 |
| N50_Rc12_L100 | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.915  | 0.9013 |   0.3486 | =     | 0.6156 |  1.5163 |
| N50_Rc12_L140 | max_CR   | NSGA-II+REPAIR | NSGA-II+RND    | 0.631  | 0.599  |   0.0159 | +     | 0.7156 |  5.3422 |

### Aggregation regime, NSGA-II+SA: median largest coverage (%) with in-degree at most d

| scenario      |   0 |   1 |   2 |   3 |   4 |   5 |
|:--------------|----:|----:|----:|----:|----:|----:|
| N30_Rc8_L100  |   7 |  37 |  49 |  49 |  49 |  49 |
| N30_Rc8_L140  |   3 |  24 |  30 |  30 |  30 |  30 |
| N30_Rc12_L100 |  12 |  56 |  71 |  71 |  71 |  71 |
| N30_Rc12_L140 |   6 |  33 |  45 |  45 |  45 |  45 |
| N50_Rc8_L100  |   7 |  40 |  69 |  69 |  69 |  69 |
| N50_Rc8_L140  |   3 |  22 |  44 |  44 |  44 |  44 |
| N50_Rc12_L100 |  12 |  57 |  93 |  93 |  93 |  93 |
| N50_Rc12_L140 |   6 |  30 |  63 |  63 |  63 |  63 |

Lifetime ratio T_agg(d)/T_agg(0): Rc=8: 1.000 / 0.503 / 0.336 / 0.252 / 0.202 / 0.168, Rc=12: 1.000 / 0.507 / 0.340 / 0.255 / 0.205 / 0.171

### Cross-regime HV loss (% of the HV of the front optimised for the evaluating regime); Holm over 16 tests per column

| scenario      | algorithm   |   agg_loss_pct |   p_agg |   raw_loss_pct |   p_raw |   p_agg_holm |   p_raw_holm |
|:--------------|:------------|---------------:|--------:|---------------:|--------:|-------------:|-------------:|
| N30_Rc12_L100 | NSGA-II     |        -0.9556 |  0.4771 |        32.4202 |       0 |       1      |            0 |
| N30_Rc12_L100 | NSGA-II+SA  |         9.4446 |  0      |        30.2094 |       0 |       0      |            0 |
| N30_Rc12_L140 | NSGA-II     |         0.0781 |  0.9354 |        30.0289 |       0 |       1      |            0 |
| N30_Rc12_L140 | NSGA-II+SA  |        10.2421 |  0      |        31.6175 |       0 |       0      |            0 |
| N30_Rc8_L100  | NSGA-II     |         2.8751 |  0.6266 |        25.8318 |       0 |       1      |            0 |
| N30_Rc8_L100  | NSGA-II+SA  |         7.8978 |  0      |        30.7924 |       0 |       0.0003 |            0 |
| N30_Rc8_L140  | NSGA-II     |         0.0017 |  0.9909 |        20.8723 |       0 |       1      |            0 |
| N30_Rc8_L140  | NSGA-II+SA  |         6.758  |  0.0003 |        32.2658 |       0 |       0.0022 |            0 |
| N50_Rc12_L100 | NSGA-II     |         7.0347 |  0      |        35.4249 |       0 |       0      |            0 |
| N50_Rc12_L100 | NSGA-II+SA  |         9.737  |  0      |        36.3673 |       0 |       0      |            0 |
| N50_Rc12_L140 | NSGA-II     |         1.0412 |  0.9354 |        36.4075 |       0 |       1      |            0 |
| N50_Rc12_L140 | NSGA-II+SA  |         8.5704 |  0      |        38.7725 |       0 |       0.0001 |            0 |
| N50_Rc8_L100  | NSGA-II     |        -0.3632 |  0.8712 |        32.6155 |       0 |       1      |            0 |
| N50_Rc8_L100  | NSGA-II+SA  |        10.5657 |  0      |        35.7339 |       0 |       0      |            0 |
| N50_Rc8_L140  | NSGA-II     |         4.4157 |  0.4123 |        27.622  |       0 |       1      |            0 |
| N50_Rc8_L140  | NSGA-II+SA  |        10.6923 |  0      |        37.456  |       0 |       0      |            0 |

Validation (aggregation): 16800 front points; T_model == T_exact (max-flow) in 100.0 %; true-distance FND >= model in 100.0 %; median ratio FND_true/T_model 1.001

Run time (s), mean per algorithm: {'AGE-MOEA-II': 36.6, 'NSGA-II': 3.8, 'NSGA-II+REPAIR': 4.7, 'NSGA-II+RND': 4.7, 'NSGA-II+SA': 4.8, 'NSGA-III': 7.9, 'SMS-EMOA': 9.3}

Budget check: n_evals == budget in 100 % of runs

### Largest coverage (median %) with lifetime >= about T_max/3 (c <= 2 raw, d <= 2 aggregation), NSGA-II+SA

| scenario      |   ratio_raw_c2 |   ratio_agg_d2 |   median_CR_raw |   median_CR_agg |
|:--------------|---------------:|---------------:|----------------:|----------------:|
| N30_Rc8_L100  |          0.335 |          0.336 |            15   |            49   |
| N30_Rc8_L140  |          0.335 |          0.336 |             7   |            30   |
| N30_Rc12_L100 |          0.336 |          0.34  |            23   |            71   |
| N30_Rc12_L140 |          0.336 |          0.34  |            13   |            45   |
| N50_Rc8_L100  |          0.335 |          0.336 |            14   |            69   |
| N50_Rc8_L140  |          0.335 |          0.336 |             7.5 |            44.5 |
| N50_Rc12_L100 |          0.336 |          0.34  |            22.5 |            93   |
| N50_Rc12_L140 |          0.336 |          0.34  |            12   |            63   |
