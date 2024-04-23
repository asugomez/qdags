#!/bin/bash
chmod a+x *.sh
# ./runqueries-$file-bfs-sorted.sh type_fun size_queue > ../../../outputs/partial/dfuds/backtracking/$file.txt
# run tests for each type_fun and each size_queue
for type_fun in {0..1}; do
  # echo type fun
  echo type_fun$type_fun>> ../../../outputs/partial/dfuds/nonFixedQueue/results.csv
  echo size_queue,j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2>> ../../../outputs/partial/dfuds/backtracking/results.csv
#  echo >> ../../../outputs/partial/dfuds/nonFixedQueue/results.csv
  printf $size_queue, >> ../../../outputs/partial/dfuds/nonFixedQueue/results.csv
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do

    input_file="./runqueries-$file-bfs-sorted.sh"
    output_file="./runqueries-$file-bfs-sorted-args.sh"
    # run query with these arguments: datasets type_fun size_queue
    # Add priorities, type_fun and size_queue
    # Iterate over each line of the input file
    while IFS= read -r line || [ -n "$line" ]; do
      modified_line="$line $type_fun"
      echo "$modified_line"
    done < "$input_file" > "$output_file"

    results_file="../../../outputs/partial/dfuds/nonFixedQueue/$file.txt"

    chmod +x $output_file

    $output_file >> $results_file
    # Calculate mean using awk
    mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
    printf $mean, >> ../../../outputs/partial/dfuds/nonFixedQueue/results.csv
  done
  echo "" >> ../../../outputs/partial/dfuds/nonFixedQueue/results.csv

done