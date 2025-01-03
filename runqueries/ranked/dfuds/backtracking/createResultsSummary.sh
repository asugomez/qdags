#!/bin/bash

# run tests for each type_fun and each k
for type_fun in 0; do
  chmod a+x *.sh
  data_csv="../../../outputs/ranked/dfuds/backtracking/results-f$type_fun-v2.csv"
  # echo type fun
  echo "type_fun : $type_fun"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
  # echo >> ../../../outputs/ranked/dfuds/backtracking/results.csv
  for k in 1 10 100 1000; do
    # echo k
    echo "k: $k"
    printf "$k;" >> $data_csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
        results_file="../../../outputs/ranked/dfuds/backtracking/$file-f$type_fun-k$k.txt"

        chmod +x $results_file

        # Calculate mean using awk
        mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
        printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done
done