#!/bin/bash

# ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/dfuds/backtracking/$file.txt
# run tests for each type_fun and each k
for type_fun in 0; do
  chmod a+x *.sh
  data_csv="../../../../outputs/query1000results/ranked/dfuds/backtracking/results-f$type_fun-v1000.csv"
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

      output_file="./runqueries-$file-bfs-sorted.sh"
      input_file="./runqueries-$file-bfs-sorted-args.sh"

      results_file="../../../../outputs/query1000results/ranked/dfuds/backtracking/$file-f$type_fun-k$k-v1000.txt"

      # Create the modified script with the updated last argument
      while IFS= read -r line || [ -n "$line" ]; do
        if [[ "$line" =~ [[:space:]]$ ]]; then
          line="${line% }"  # Remove trailing space if present
        fi
        #modified_line=$(echo "$line" | sed -E "s/[0-9]+$/ $k/") # Replace last number with $k
        #modified_line=$(echo "$line" | sed -E "s/[[:space:]]([[:digit:]]+)$/\1$k/")
        # Correctly replace the last argument and remove double spaces
        modified_line=$(echo "$line" | sed -E "s/[[:space:]]([[:digit:]]+)$/ $k/" | sed -E "s/  +/ /g")
        echo "$modified_line"
      done < "$input_file" > "$output_file"

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
