#!/bin/bash

cd "$(dirname "$0")"

echo "================================"
echo "Running All 9 Controlled Experiments"
echo "τ = 0.10 (frozen)"
echo "Wall-clock budget: 60 seconds each"
echo "Total estimated time: ~10 minutes"
echo "================================"
echo ""

START_TIME=$(date +%s)

for policy in rrt quadtree_no_energy quadtree_energy; do
  for seed in 42 43 44; do
    echo "========================================="
    echo "Running: ${policy} seed=${seed}"
    echo "========================================="
    ./build/src/experiment_runner configs/${policy}_seed${seed}.cfg
    echo ""
  done
done

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo "================================"
echo "All experiments complete!"
echo "Total time: ${ELAPSED} seconds"
echo "Results saved in results/ directory"
echo "================================"
