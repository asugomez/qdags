#!/bin/bash
chmod a+x *.sh
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    # get the number of datasets for each query
    echo "file : $file"
    input_file="./runqueries-$file-bfs-sorted.sh"
    count=1
    while IFS= read -r line || [ -n "$line" ]; do
      read -r t1 t2 t3 t4 t5 <<< $line
#      echo "t1: $t1" # ./j4
#      echo "t2: $t2"
#      echo "t3: $t3"
#      echo "t4: $t4"
#      echo "t5: $t5"

      ../../data/priorities/createRandomPriorities.sh $t2 "$file/pri1-$count"
      ../../data/priorities/createRandomPriorities.sh $t3 "$file/pri2-$count"
      ../../data/priorities/createRandomPriorities.sh $t4 "$file/pri3-$count"
      ../../data/priorities/createRandomPriorities.sh $t5 "$file/pri4-$count"
      count=$(($count + 1))
#      echo "count : $count"
    done < $input_file

done
