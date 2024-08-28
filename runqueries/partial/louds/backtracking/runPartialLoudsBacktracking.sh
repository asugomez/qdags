#!/bin/bash

# ./runqueries-$file-bfs-sorted.sh type_fun k > ../../../outputs/partial/louds/backtracking/$file.txt
# run tests for each type_fun and each k
for type_fun in 0; do
  chmod a+x *.sh
  data_csv="../../../outputs/partial/louds/backtracking/results-f$type_fun.csv"
  # echo type fun
  echo "type_fun : $type_fun"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
  # echo >> ../../../outputs/partial/louds/backtracking/results.csv
  for k in 1 10 100 1000; do
    # echo k
    echo "k: $k"
    printf "$k;" >> $data_csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      echo "file: $file"
      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"
      # run query with these arguments: datasets type_fun k
      # Add type_fun and k
      # Iterate over each line of the input file
      while IFS= read -r line || [ -n "$line" ]; do
        if [[ "$line" =~ [[:space:]]$ ]]; then
          line="${line% }"  # Remove the trailing space
        fi
        modified_line="$line $type_fun $k"
        echo "$modified_line"
      done < "$input_file" > "$output_file"

      results_file="../../../outputs/partial/louds/backtracking/$file-f$type_fun-k$k.txt"

      chmod +x $output_file

      $output_file >> $results_file
      # Calculate mean using awk
      mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
      printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done
done