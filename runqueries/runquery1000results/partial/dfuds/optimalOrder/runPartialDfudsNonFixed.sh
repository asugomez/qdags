#!/bin/bash

# ./runqueries-$file-bfs-sorted.sh type_fun k > ../../../outputs/partial/dfuds/optimalOrder/$file.txt
# run tests for each type_fun and each k
for type_fun in 0; do
  chmod a+x *.sh
  data_csv="../../../../outputs/query1000results/partial/dfuds/optimalOrder/results-f$type_fun-v1000-time.csv"
  # echo type fun
  echo "type_fun : $type_fun"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" >> $data_csv
  for k in 1 10 100 1000; do
    # echo k
    echo "k: $k"
    printf "$k;" >> $data_csv
    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      echo "file: $file"
      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"
      # run query with these arguments: datasets type_fun k
      # Add ype_fun and k
      # Iterate over each line of the input file
      while IFS= read -r line || [ -n "$line" ]; do
        if [[ "$line" =~ [[:space:]]$ ]]; then
          line="${line% }"  # Remove the trailing space
        fi
        modified_line="$line $type_fun $k"
        echo "$modified_line"
      done < "$input_file" > "$output_file"

      results_file="../../../../outputs/query1000results/partial/dfuds/optimalOrder/$file-f$type_fun-k$k-v1000-time.txt"

      chmod +x $output_file

      $output_file >> $results_file
      # Calculate the mean of hexadecimal numbers
      sum=0
      count=0

      while IFS= read -r hex_number || [ -n "$hex_number" ]; do
        # Convert the hexadecimal number to decimal
        decimal_value=$((0x$hex_number))
        sum=$((sum + decimal_value))
        count=$((count + 1))
      done < "$results_file"

      # Calculate mean and handle division by zero
      if [ $count -gt 0 ]; then
        mean=$((sum / count))
      else
        mean=0
      fi

      printf "$mean;" >> $data_csv
    done
    echo "" >> $data_csv
  done

done