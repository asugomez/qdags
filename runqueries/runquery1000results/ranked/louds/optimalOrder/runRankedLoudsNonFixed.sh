#!/bin/bash
# run tests for each type_fun and each k
for type_fun in 0 1; do
  chmod a+x *.sh
  time_csv="../../../../outputs/query1000results/ranked/louds/optimalOrder/results-f$type_fun-time.csv"
  nodes_csv="../../../../outputs/query1000results/ranked/louds/optimalOrder/results-f$type_fun-nodes.csv"

  echo "type_fun : $type_fun"

  # CSV headers
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" > "$time_csv"
  echo "k;j3;j4;p2;p3;p4;s1;s2;s3;s4;t2;t3;t4;ti2;ti3;ti4;tr1;tr2" > "$nodes_csv"

  for k in 1 10 100 1000; do
    # echo k
    echo "k: $k"
    printf "$k;" >> "$time_csv"
    printf "$k;" >> "$nodes_csv"

    for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
      # get the number of datasets for each query
      echo "file: $file"

      input_file="./runqueries-$file-bfs-sorted.sh"
      output_file="./runqueries-$file-bfs-sorted-args.sh"


      # Create the modified script with the updated last argument
      while IFS= read -r line || [ -n "$line" ]; do
        modified_line=$(awk -v k="$k" '{ $NF = k; print }' <<< "$line")
        echo "$modified_line"
      done < "$input_file" > "$output_file"

      results_file="../../../../outputs/query1000results/ranked/louds/optimalOrder/$file-f$type_fun-k$k-v1000-results.txt"
      chmod +x $output_file
      > "$results_file"

      $output_file >> $results_file

      # Crear archivos temporales para nodos (a1 = hex, odd lines) y tiempo (a2 = float, even lines)
      nodes_file=$(mktemp)
      times_file=$(mktemp)

      awk 'NR % 2 == 1 { print }' "$results_file" > "$nodes_file"
      awk 'NR % 2 == 0 { print }' "$results_file" > "$times_file"

      ## Calcular media de nodos (hexadecimal)
      sum=0
      count=0
      while IFS= read -r hex_number || [ -n "$hex_number" ]; do
        decimal_value=$((16#$hex_number))
        sum=$((sum + decimal_value))
        count=$((count + 1))
      done < "$nodes_file"

      if [ $count -gt 0 ]; then
        mean_nodes=$((sum / count))
      else
        mean_nodes=0
      fi

      ## Calcular media de tiempo (flotante)
      mean_time=$(awk '{sum += $1} END {if (NR > 0) printf "%.6f", sum / NR; else print 0}' "$times_file")

      ## Escribir a CSV
      printf "$mean_time;" >> "$time_csv"
      printf "$mean_nodes;" >> "$nodes_csv"

      rm "$nodes_file" "$times_file"
    done

    echo "" >> "$time_csv"
    echo "" >> "$nodes_csv"
  done
done


