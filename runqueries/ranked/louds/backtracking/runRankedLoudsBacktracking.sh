#!/bin/bash
# ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/louds/backtracking/$file.txt
# run tests for each type_fun and each size_queue
for((type_fun = 0; type_fun < 2; type_fun +=1)); do
  # echo type fun
  echo type_fun$type_fun >> ../../../outputs/ranked/louds/backtracking/results.csv
  echo size_queue,j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2>> ../../../outputs/ranked/louds/backtracking/results.csv
#  echo >> ../../../outputs/ranked/louds/backtracking/results.csv
  for((size_queue = 10; size_queue <= 1010; size_queue += 100)); do
    # echo size_queue
    printf $size_queue, >> ../../../outputs/ranked/louds/backtracking/results.csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      # Obtener los argumentos del runqueries
      # obtengo los argumentos de la linea 10 del archivo runqueries-$file-bfs-sorted.sh
      linea=$(awk 'NR==10 { print $2, $3, $4, $5 }' ./runqueries-$file-bfs-sorted.sh)
      # Asignar los valores a variables separadas
      read t1 t2 t3 t4 <<< "$linea"
      echo t1 $t1 t2 $t2 t3 $t3 t4 $t4
      # Create priorities

      priority_file_1="../../../../data/priorities/p1"
      priority_file_2="../../../../data/priorities/p2"
      priority_file_3="../../../../data/priorities/p3"
      priority_file_4="../../../../data/priorities/p4"

      ../../../../data/priorities/createRandomPriorities.sh $t1 1 > $priority_file_1
      ../../../../data/priorities/createRandomPriorities.sh $t2 2 > $priority_file_2
      ../../../../data/priorities/createRandomPriorities.sh $t3 3 > $priority_file_3
      ../../../../data/priorities/createRandomPriorities.sh $t4 4 > $priority_file_4

      echo name $priority_file_1 $priority_file_2 $priority_file_3 $priority_file_4

      archivo="../../../outputs/ranked/louds/backtracking/$file.txt"
      # run query with these arguments: datasets type_fun size_queue
      # Check if the i-th argument is emtpy
      if [ -z "$t2" ]; then
        ./runqueries-$file-bfs-sorted.sh priority_file_1 type_fun size_queue > $archivo
      elif [ -z "$t3" ]; then
        ./runqueries-$file-bfs-sorted.sh priority_file_1 priority_file_2 type_fun size_queue > $archivo
      elif [ -z "$t4" ]; then
        ./runqueries-$file-bfs-sorted.sh priority_file_1 priority_file_2 priority_file_3 type_fun size_queue > $archivo
      else
        ./runqueries-$file-bfs-sorted.sh priority_file_1 priority_file_2 priority_file_3 priority_file_4  type_fun size_queue > $archivo
      fi
      # Calcular el promedio utilizando awk
      promedio=$(awk '{ suma += $1 } END { print suma / NR }' "$archivo")
      echo promedio $promedio
      printf $promedio, >> ../../../outputs/ranked/louds/backtracking/results.csv
    done
    echo  >> ../../../outputs/ranked/louds/backtracking/results.csv
  done
done