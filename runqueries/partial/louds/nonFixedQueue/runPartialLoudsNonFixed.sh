#!/bin/bash

# ./runqueries-$file-bfs-sorted.sh type_fun > ../../../outputs/partial/louds/nonFixedQueue/$file.txt
# run tests for each type_fun and each size_queue
for((type_fun = 0; type_fun < 2; type_fun +=1)); do
  # echo type fun
  echo type_fun$type_fun >> ../../../outputs/partial/louds/nonFixedQueue/results.csv
  echo j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2>> ../../../outputs/partial/louds/nonFixedQueue/results.csv
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    # run query with these arguments: datasets type_fun size_queue
    archivo="../../../outputs/partial/louds/nonFixedQueue/$file.txt"
    ./runqueries-$file-bfs-sorted.sh size_queue > $archivo
    # Calcular el promedio utilizando awk
    promedio=$(awk '{ suma += $1 } END { print suma / NR }' "$archivo")
    echo promedio $promedio
    printf $promedio, >> ../../../outputs/partial/louds/nonFixedQueue/results.csv
  done
  echo  >> ../../../outputs/partial/louds/nonFixedQueue/results.csv
done