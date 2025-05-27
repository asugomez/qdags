#!/bin/bash

# Define directories relative to ./ranked
base_dirs=(
  "../../outputs/query1000results/ranked/dfuds/backtracking"
  "../../outputs/query1000results/ranked/dfuds/optimalOrder"
  "../../outputs/query1000results/ranked/louds/backtracking"
  "../../outputs/query1000results/ranked/louds/optimalOrder"
)

# Define files and k values
files=(j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2)
k_values=(1)

for dir in "${base_dirs[@]}"; do
  for k in "${k_values[@]}"; do
    for file in "${files[@]}"; do
      results_file="$dir/${file}-f0-k${k}-v1000-results.txt"
      nodes_file="$dir/${file}-f0-k${k}-v1000-nodes.txt"
      time_file="$dir/${file}-f0-k${k}-v1000-time.txt"

      if [ -f "$results_file" ]; then
        awk 'NR % 2 == 1' "$results_file" > "$nodes_file"
        awk 'NR % 2 == 0' "$results_file" > "$time_file"
        echo "✓ Processed: $results_file → nodes & time"
      else
        echo "✗ Missing: $results_file"
      fi
    done
  done
done
