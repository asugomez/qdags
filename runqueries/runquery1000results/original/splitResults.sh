#!/bin/bash

for k in 1; do
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    #results_file="../../outputs/query1000results/original/${file}-k${k}-results-gradual.txt"
    results_file="../../outputs/query1000results/original/${file}-results-ranked.txt"
    nodes_file="../../outputs/query1000results/original/${file}-nodes-ranked.txt"
    time_file="../../outputs/query1000results/original/${file}-time-ranked.txt"

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
