#!/bin/bash
chmod a+x *.sh
data_csv="../outputs/original/results-time.csv"

echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
for k in 10000 100000 1000000; do # 1 10 100 1000
  # echo k
  echo "k: $k"
  printf "$k;" >> $data_csv
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    echo "file: $file"
    input_file="./runqueries-$file-bfs-sorted.sh"
    output_file="./runqueries-$file-bfs-sorted-args.sh"
    results_file="../outputs/original/$file-k$k-time.txt"

    while IFS= read -r line || [ -n "$line" ]; do
      if [[ "$line" =~ [[:space:]]$ ]]; then
        line="${line% }"  # Remove the trailing space
      fi
      modified_line="$line $k"
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

    # Calculate mean using awk
    #mean=$(awk '{ suma += $(('0'x1)) } END { print suma / NR }' "../outputs/original/test")
    printf "$mean;" >> $data_csv
  done
  echo "" >> $data_csv
done

