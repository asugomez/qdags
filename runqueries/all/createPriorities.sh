#!/bin/bash
chmod a+x *.sh
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
    # get the number of datasets for each query
    read -r t1 t2 t3 t4 <<< "$(awk 'NR==1 {print $2 " " $3 " " $4 " " $5 }' ./runqueries-$file-bfs-sorted.sh)"

    echo "t1: $t1"
    echo "t2: $t2"
    echo "t3: $t3"
    echo "t4: $t4"
    # Create priorities
    priority_file_1="../../data/priorities/pri1-$file"
    priority_file_2="../../data/priorities/pri2-$file"
    priority_file_3="../../data/priorities/pri3-$file"
    priority_file_4="../../data/priorities/pri4-$file"

    ../../data/priorities/createRandomPriorities.sh $t1 "pri1-$file"
    ../../data/priorities/createRandomPriorities.sh $t2 "pri2-$file"
    ../../data/priorities/createRandomPriorities.sh $t3 "pri3-$file"
    ../../data/priorities/createRandomPriorities.sh $t4 "pri4-$file"

done
