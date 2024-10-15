#!/bin/bash

process_files() {
  local directory=$1

  input_file="$directory/results-f0.csv"

  chmod +x $input_file

  index=1
  for k in 1 10 100 1000; do
    results_file="$directory/s1-f0-k$k.txt"

    echo "k: $k"

    # Calculate mean using awk
    mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")

    echo "mean: $mean"

    # Use awk to modify the second row and seventh column
    awk -v new_value="$mean" -v row="$index" 'BEGIN { FS=OFS="," }
    NR==row {$7=new_value} 1' "$input_file" > tmpfile && mv tmpfile "$FILE"

    index=$((index + 1))  # Increment the index
  done
}

process_files "$1"