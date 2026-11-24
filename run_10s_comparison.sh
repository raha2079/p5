#!/bin/bash
cd /Users/gray/Documents/AI\ robotics/project\ 5

echo "Running 10s comparison experiments (9 runs)..."
for seed in 42 43 44; do
  for planner in quadtree_energy quadtree_no_energy rrt; do
    echo "========================================="
    echo "Running: ${planner} seed ${seed}"
    ./build/src/experiment_runner configs/comp_10s_${planner}_seed${seed}.cfg 2>&1 | grep -E "(Planner|Seed|States:|Iterations:|Attempted edges:|Accepted edges:|Acceptance rate:|Energy -|Leaf cells:|Efficiency:)" | head -15
    echo ""
  done
done
echo "All 10s comparison experiments complete!"
