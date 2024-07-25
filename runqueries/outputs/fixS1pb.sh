#!/bin/bash

process_files() {
  local directory=$1

  for k in 1 10 100 1000; do
    input_file="$directory/s1-f0-k$k.txt"
    output_file="$directory/s1-f0-k$k-fixed.txt"

    chmod +x "$input_file"
    chmod +x "$output_file"

    awk '{print $2}' "$input_file" > "$output_file"
    mv "$output_file" "$input_file"
  done
}

# Call the function with the directory name as an argument
process_files "$1"
