#!/bin/sh
# Published repair baselines (review round 4, item 7). Run from the repository root after run_extra.sh.
set -e
A1_CONFIG=config_pub.json A1_RES=results_pub python3 experiments/run_a1.py --exp PB --workers 4
A1_CONFIG=config_pub_agg.json A1_RES=results_pub_agg python3 experiments/run_a1.py --exp PA --workers 4
echo EXTRA2_DONE
