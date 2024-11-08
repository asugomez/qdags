#!/bin/bash
chmod a+x *.sh
data_csv="../../outputs/query1000results/lqdags/results-v1000.csv"

echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
for k in 1 10 100 1000; do
  # echo k
  echo "k: $k"
  printf "$k;" >> $data_csv
  for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    echo "file: $file"
    input_file="./runqueries-$file-bfs-sorted.sh"
    output_file="./runqueries-$file-bfs-sorted-args.sh"
    results_file="../../outputs/query1000results/lqdags/$file-k$k-v1000.txt"

    while IFS= read -r line || [ -n "$line" ]; do
      if [[ "$line" =~ [[:space:]]$ ]]; then
        line="${line% }"  # Remove the trailing space
      fi
      modified_line="$line $k"
      echo "$modified_line"
    done < "$input_file" > "$output_file"

    chmod +x $output_file

    $output_file >> $results_file
    # Calculate mean using awk
    mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
    printf "$mean;" >> $data_csv
  done
  echo "" >> $data_csv
done

