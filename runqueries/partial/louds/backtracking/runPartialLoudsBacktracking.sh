#!/bin/bash
chmod a+x *.sh
data_csv="../../../outputs/partial/louds/backtracking/results.csv"
echo "" > $data_csv

# ./runqueries-$file-bfs-sorted.sh type_fun size_queue > ../../../outputs/partial/louds/backtracking/$file.txt
# run tests for each type_fun and each size_queue
for type_fun in {0..1}; do
  # echo type fun
  echo "type_fun;$type_fun" >> $data_csv
  echo "size_queue;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
#  echo >> ../../../outputs/partial/louds/backtracking/results.csv
  for size_queue in 1 10 ; do #100 1000; do
    # echo size_queue
    printf "$size_queue;" >> $data_csv
    for file in j3; do #j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do

      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"
      # run query with these arguments: datasets type_fun size_queue
      # Add priorities; type_fun and size_queue
      # Iterate over each line of the input file
      while IFS= read -r line || [ -n "$line" ]; do
        modified_line="$line $type_fun $size_queue"
        echo "$modified_line"
      done < "$input_file" > "$output_file"

      results_file="../../../outputs/partial/louds/backtracking/$file-f$type_fun-s$size_queue.txt"

      chmod +x $output_file

      $output_file >> $results_file
      # Calculate mean using awk
      mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
      printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done
done