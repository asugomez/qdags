#!/bin/bash
chmod a+x *.sh
# ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/louds/backtracking/$file.txt
# run tests for each type_fun and each size_queue
for((type_fun = 0; type_fun < 2; type_fun +=1)); do
  # echo type fun
  echo type_fun,$type_fun >> ../../../outputs/ranked/louds/backtracking/results.csv
  echo size_queue,j3,j4,p2,p3,p4,s1,s2,s3,s4,t2,t3,t4,ti2,ti3,ti4,tr1,tr2>> ../../../outputs/ranked/louds/backtracking/results.csv
  for size_queue in 1 10 100 1000; do
    # echo size_queue
    printf $size_queue, >> ../../../outputs/ranked/louds/backtracking/results.csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      # get the number of datasets for each query
      line=$(awk 'NR==10 { print $2, $3, $4, $5 }' ./runqueries-$file-bfs-sorted.sh)
      read t1 t2 t3 t4 <<< "$line"

      # Create priorities
      priority_file_1="../../../../data/priorities/pri1"
      priority_file_2="../../../../data/priorities/pri2"
      priority_file_3="../../../../data/priorities/pri3"
      priority_file_4="../../../../data/priorities/pri4"

      dataset1="../../../../data/all/"$(../../../../data/chooseDataset.sh)
      dataset2="../../../../data/all/"$(../../../../data/chooseDataset.sh)
      dataset3="../../../../data/all/"$(../../../../data/chooseDataset.sh)
      dataset4="../../../../data/all/"$(../../../../data/chooseDataset.sh)

      ../../../../data/priorities/createRandomPriorities.sh $dataset1 "pri1"
      ../../../../data/priorities/createRandomPriorities.sh $dataset2 "pri2"
      ../../../../data/priorities/createRandomPriorities.sh $dataset3 "pri3"
      ../../../../data/priorities/createRandomPriorities.sh $dataset4 "pri4"

      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"

      # Add priorities, type_fun and size_queue
      # Iterate over each line of the input file
      while IFS= read -r line || [ -n "$line" ]; do
        # Append priorities, type_fun and size_queue to the end of the line
        modified_line=""
        # Check if the i-th argument is emtpy
        if [ -z "$t2" ]; then # 1 dataset
          modified_line="$line $priority_file_1 $type_fun $size_queue"
        elif [ -z "$t3" ]; then # 2 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $type_fun $size_queue"
        elif [ -z "$t4" ]; then # 3 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $priority_file_3 $type_fun $size_queue"
        else # 4 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $priority_file_3 $priority_file_4 $type_fun $size_queue"
        fi
        echo "$modified_line"
      done < "$input_file" > "$output_file"

      results_file="../../../outputs/ranked/louds/backtracking/$file-f$type_fun-s$size_queue.txt"

      chmod +x $output_file

      $output_file >> $results_file

      # Calculate mean using awk
      mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
      printf $mean, >> ../../../outputs/ranked/louds/backtracking/results.csv
    done
    echo "" >> ../../../outputs/ranked/louds/backtracking/results.csv
  done
done