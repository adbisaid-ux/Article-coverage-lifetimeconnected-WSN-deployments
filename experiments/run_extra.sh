#!/bin/sh
# Extra campaigns requested by review round 4 (items 10 and 16). Run from the repository root.
set -e
A1_CONFIG=config_var.json A1_RES=results_var python3 experiments/run_a1.py --exp V --workers 4
A1_CONFIG=config_tune.json A1_RES=results_tune python3 experiments/run_a1.py --exp TU --workers 4
echo EXTRA_DONE
