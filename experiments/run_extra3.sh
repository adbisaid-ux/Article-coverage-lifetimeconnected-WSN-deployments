#!/bin/sh
# Rate sensitivity of the repairs (review round 5, item 6). Run from the repository root.
set -e
A1_CONFIG=config_pubrate.json A1_RES=results_pubrate python3 experiments/run_a1.py --exp PR --workers 4
echo EXTRA3_DONE
