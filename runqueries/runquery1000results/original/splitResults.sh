#!/bin/bash

for k in 1; do
  for file in t4 ti4; do
    results_file="../../outputs/query1000results/original/${file}-k${k}-results-gradual.txt"
    nodes_file="../../outputs/query1000results/original/${file}-k${k}-nodes.txt"
    time_file="../../outputs/query1000results/original/${file}-k${k}-time.txt"

    if [ -f "$results_file" ]; then
      # Extract odd lines (1, 3, 5, ...) → nodes
      awk 'NR % 2 == 1' "$results_file" > "$nodes_file"
      # Extract even lines (2, 4, 6, ...) → time
      awk 'NR % 2 == 0' "$results_file" > "$time_file"
      echo "Processed $results_file → $nodes_file & $time_file"
    else
      echo "File not found: $results_file"
    fi
  done
done
