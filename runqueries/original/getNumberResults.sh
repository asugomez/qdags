#!/bin/bash
chmod a+x *.sh

results_number_csv="../outputs/original/results-number.csv"

# CSV header
echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" > "$results_number_csv"

for k in 1; do
  echo "k: $k"
  printf "$k;" >> "$results_number_csv"

  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    echo "file: $file"
    input_file="./runqueries-$file-bfs-sorted.sh"
    output_file="./runqueries-$file-bfs-sorted-args.sh"

    # Add argument k to each line
    while IFS= read -r line || [ -n "$line" ]; do
      line="${line%" "}"
      echo "$line $k"
    done < "$input_file" > "$output_file"

    temp_file=$(mktemp)
    chmod +x "$output_file"
    "$output_file" > "$temp_file"

    # Keep only odd lines in results_file
    results_file="../outputs/original/$file-results-number.txt"
    awk 'NR % 2 == 1' "$temp_file" > "$results_file"
    rm "$temp_file"

    ## Compute average number of nodes (hex from odd lines)
    sum=0
    count=0
    while IFS= read -r hex_number; do
      decimal_value=$((16#$hex_number))
      sum=$((sum + decimal_value))
      count=$((count + 1))
    done < "$results_file"

    if [ $count -gt 0 ]; then
      mean_nodes=$((sum / count))
    else
      mean_nodes=0
    fi

    ## Write to CSV
    printf "$mean_nodes;" >> "$results_number_csv"
  done

  echo "" >> "$results_number_csv"
done
