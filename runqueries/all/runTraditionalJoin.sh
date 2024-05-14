#!/bin/bash
chmod a+x *.sh
data_csv="../outputs/all/results.csv"
echo "" > $data_csv

echo "j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2">> $data_csv
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  echo "file: $file"
  input_file="./runqueries-$file-bfs-sorted.sh"

  results_file="../outputs/all/$file.txt"

  chmod +x $input_file

  $input_file >> $results_file
  # Calculate mean using awk
  mean=$(awk '{ suma += $1 } END { print suma / NR }' "$results_file")
  printf "$mean;" >> $data_csv
done
echo "" >> $data_csv

