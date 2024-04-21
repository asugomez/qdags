#!/bin/bash
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  ./runqueries-$file-bfs-sorted-fun-0.sh > ../../../outputs/partial/louds/backtracking/$file.txt
  ./runqueries-$file-bfs-sorted-fun-1.sh > ../../../outputs/partial/louds/backtracking/$file.txt
done