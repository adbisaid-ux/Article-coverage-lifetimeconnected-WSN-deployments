# Coverage versus lifetime in connected wireless sensor deployments

Code, configurations, seeds, raw results and analysis scripts for the manuscript
"Coverage versus lifetime in connected wireless sensor deployments: a load-balanced lifetime model and a
connectivity-preserving operator for multi-objective optimisation" (S. Adbi, H. Mouncif).

Python 3.13; exact library versions in `requirements.txt` (also recorded in each `results*/environment.json`).

## Layout
- `src/model.py`: problem model (connected coverage; lifetime without aggregation, `c_tree`, and with
  aggregation, `agg_tree_assign`), exact references (`exact_c`, `exact_agg_degree`) and the operator `agent_step`
  (SA-op; mode 1 = random moves, mode 2 = repair only).
- `src/repairs_pub.py`: the published repairs KY-RH (Konstantinidis and Yang, 2011) and DD-MCSA (Deghbouch and
  Debbat, 2022), re-implemented as baselines.
- `src/moea.py`: NSGA-II(+SA), DE baselines (DE-CR, EPS-DE), pymoo wrappers, hypervolume.
- `experiments/run_a1.py`: runs a campaign (resumable). `A1_CONFIG` selects the configuration file and `A1_RES`
  the output folder; `--exp X` selects `scenarios_X` / `algorithms_X` of that configuration.
- `experiments/config*.json`: scenarios, budgets, seeds and algorithm lists of every experiment.
- `tests/`: connectivity of SA-op and optimality of the aggregation tree (`python3 tests/<file>.py`).

## Reproducing the experiments
| Experiment | Command | Output |
|---|---|---|
| Main campaign and ablation | `python3 experiments/run_a1.py --workers 4` then `--ablation` | `results/` |
| E1-E3 (controls, generality, robustness) | `A1_CONFIG=config_boucle.json A1_RES=results_boucle python3 experiments/run_a1.py --exp E12` and `--exp E3` | `results_boucle/` |
| R1, R2 (S), R3 | `A1_CONFIG=config_rev1.json A1_RES=results_rev1 python3 experiments/run_a1.py --exp R1` (then `S`, `R3`) | `results_rev1/` |
| L (5e5 evaluations) | `A1_CONFIG=config_long.json A1_RES=results_long python3 experiments/run_a1.py --exp L` | `results_long/` |
| AG, XR (aggregation regime) | `A1_CONFIG=config_agg.json A1_RES=results_agg python3 experiments/run_a1.py --exp AG` (then `XR`) | `results_agg/` |
| V (variability across algorithm seeds) | `A1_CONFIG=config_var.json A1_RES=results_var python3 experiments/run_a1.py --exp V` | `results_var/` |
| TU, TN (tuned NSGA-II) | `A1_CONFIG=config_tune.json A1_RES=results_tune python3 experiments/run_a1.py --exp TU`, then `A1_CONFIG=config_tuned.json A1_RES=results_tuned ... --exp TN` | `results_tune/`, `results_tuned/` |
| PB, PA (published repairs KY-RH and DD-MCSA) | `A1_CONFIG=config_pub.json A1_RES=results_pub python3 experiments/run_a1.py --exp PB`; `A1_CONFIG=config_pub_agg.json A1_RES=results_pub_agg ... --exp PA` | `results_pub/`, `results_pub_agg/` |
| PR (repair rates 0.5 and 1) and repair displacement | `sh experiments/run_extra3.sh`; `python3 experiments/repair_moves.py` | `results_pubrate/` |
| Sensing energy (re-scoring, no new runs) | `python3 experiments/analyze_sensing.py` | `results_sensing/` |

Runs are deterministic for given seeds: re-running a configuration reproduces the stored hypervolumes exactly.

## Tables and figures
`experiments/analyze_a1.py`, `tradeoff_a1.py`, `analyze_boucle.py`, `analyze_rev1_existing.py`,
`simulate_rounds_rev1.py`, `analyze_rev1_runs.py`, `analyze_agg.py`, `analyze_extra.py` and `analyze_sensing.py` read the raw CSV files and write the tables
(`results*/tables*/`) and figures (`figures*/`) used in the paper.

## Correction history
`results_rev1/superseded_attach_0.8999/` and `results_long/superseded_attach_0.8999/` keep, for the record, runs
made with a rounding error in the attachment radius (see the README there). They are not used in the paper.

## Licence
Code: MIT licence (see `LICENSE`). Data: CC BY 4.0.
