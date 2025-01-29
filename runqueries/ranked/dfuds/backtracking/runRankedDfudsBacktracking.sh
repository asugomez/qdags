#!/bin/bash

# ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/dfuds/backtracking/$file.txt
# run tests for each type_fun and each k
for type_fun in 1; do
  chmod a+x *.sh
  data_csv="../../../outputs/ranked/dfuds/backtracking/results-f$type_fun-time.csv"
  # echo type fun
  echo "type_fun : $type_fun"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" >> $data_csv
  for k in 1 10 100 1000; do
    # echo k
    echo "k: $k"
    printf "$k;" >> $data_csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      # get the number of datasets for each query
      echo "file: $file"
      read -r t1 t2 t3 t4 <<< "$(awk 'NR==1 {print $2 " " $3 " " $4 " " $5 }' ./runqueries-$file-bfs-sorted.sh)"

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

      results_file="../../../outputs/ranked/dfuds/backtracking/$file-f$type_fun-k$k-time.txt"

      chmod +x $output_file

      $output_file >> $results_file

      # Calculate the mean of hexadecimal numbers
#      sum=0
#      count=0
#
#      while IFS= read -r hex_number || [ -n "$hex_number" ]; do
#        # Convert the hexadecimal number to decimal
#        decimal_value=$((0x$hex_number))
#        sum=$((sum + decimal_value))
#        count=$((count + 1))
#      done < "$results_file"
#
#      # Calculate mean and handle division by zero
#      if [ $count -gt 0 ]; then
#        mean=$((sum / count))
#      else
#        mean=0
#      fi

      # Calculate mean using awk
      mean=$(awk '{ suma += $(('0'x1)) } END { print suma / NR }' "$results_file")
      printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done
done