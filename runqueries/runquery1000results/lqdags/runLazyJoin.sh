#!/bin/bash
chmod a+x *.sh
data_csv="../../outputs/query1000results/lqdags/results-v1000-time-v2.csv"

echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
for k in 1 10 100 1000; do
  # echo k
  echo "k: $k"
  printf "$k;" >> $data_csv
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    echo "file: $file"
    input_file="./runqueries-$file-bfs-sorted.sh"
    output_file="./runqueries-$file-bfs-sorted-args.sh"
    results_file="../../outputs/query1000results/lqdags/$file-k$k-v1000-time-v2.txt"

    while IFS= read -r line || [ -n "$line" ]; do
      if [[ "$line" =~ [[:space:]]$ ]]; then
        line="${line% }"  # Remove the trailing space
      fi
      modified_line="$line $k"
      echo "$modified_line"
    done < "$input_file" > "$output_file"
#
    chmod +x $output_file

    $output_file >> $results_file


    # OPTION 1: NODES
    # Calculate the mean of hexadecimal numbers
#    sum=0
#    count=0

#    while IFS= read -r hex_number || [ -n "$hex_number" ]; do
#      decimal_value=$((0x$hex_number))  # Convert to decimal
#      if [ "$decimal_value" -ne 0 ]; then
#        sum=$((sum + decimal_value))
#        count=$((count + 1))
#      fi
#    done < "$results_file"

    # Calculate mean and handle division by zero
#    if [ $count -gt 0 ]; then
#      mean=$((sum / count))
#    else
#      echo "count is $count"
#      mean=0
#    fi

    # OPTION 2: time
    mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
    printf "$mean;" >> $data_csv
  done
  echo "" >> $data_csv
done

