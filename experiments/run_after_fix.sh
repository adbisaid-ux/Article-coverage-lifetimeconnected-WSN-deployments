#!/bin/sh
# Re-runs after the attachment-radius fix (18/10/2026), then the aggregation-regime campaign.
cd "$(dirname "$0")/.."
A1_CONFIG=config_rev1.json A1_RES=results_rev1 python3 experiments/run_a1.py --workers 4 --exp R3
A1_CONFIG=config_rev1.json A1_RES=results_rev1 python3 experiments/run_a1.py --workers 4 --exp S
A1_CONFIG=config_agg.json A1_RES=results_agg python3 experiments/run_a1.py --workers 4 --exp AG
A1_CONFIG=config_agg.json A1_RES=results_agg python3 experiments/run_a1.py --workers 4 --exp XR
A1_CONFIG=config_long.json A1_RES=results_long python3 experiments/run_a1.py --workers 4 --exp L
