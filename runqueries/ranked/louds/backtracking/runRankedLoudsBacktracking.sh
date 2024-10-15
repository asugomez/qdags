#!/bin/bash
# ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/louds/backtracking/$file.txt
# run tests for each type_fun and each k
for type_fun in 0; do
  chmod a+x *.sh
  data_csv="../../../outputs/ranked/louds/backtracking/results-f$type_fun.csv"
  # echo type fun
  echo "type_fun : $type_fun"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" >> $data_csv
  for k in 1 10 100 1000; do
    # echo k
    echo "size queue: $k"
    printf "$k;" >> $data_csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      #get the number of datasets for each query
      echo "file: $file"
      read -r l1 l2 l3 l4 <<< "$(awk 'NR==1 {print $2 " " $3 " " $4 " " $5 }' ./runqueries-$file-bfs-sorted.sh)"

      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"

      # Add priorities; type_fun and k
      # Iterate over each line of the input file
      count=1
      while IFS= read -r line || [ -n "$line" ]; do
        if [[ "$line" =~ [[:space:]]$ ]]; then
          line="${line% }"  # Remove the trailing space
        fi
        # Append priorities; type_fun and k to the end of the line
        priority_file_1="../../../../data/priorities/$file/pri1-$count"
        priority_file_2="../../../../data/priorities/$file/pri2-$count"
        priority_file_3="../../../../data/priorities/$file/pri3-$count"
        priority_file_4="../../../../data/priorities/$file/pri4-$count"
        modified_line=""
        # Check if the i-th argument is emtpy
        if [ -z "$t2" ]; then # 1 dataset
          modified_line="$line $priority_file_1 $type_fun $k"
        elif [ -z "$t3" ]; then # 2 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $type_fun $k"
        elif [ -z "$t4" ]; then # 3 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $priority_file_3 $type_fun $k"
        else # 4 datasets
          modified_line="$line $priority_file_1 $priority_file_2 $priority_file_3 $priority_file_4 $type_fun $k"
        fi
        echo "$modified_line"
        count=$(($count + 1))
      done < "$input_file" > "$output_file"

      results_file="../../../outputs/ranked/louds/backtracking/$file-f$type_fun-k$k.txt"

      chmod +x $output_file

      $output_file >> $results_file

      # Calculate mean using awk
      mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
      printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done
done