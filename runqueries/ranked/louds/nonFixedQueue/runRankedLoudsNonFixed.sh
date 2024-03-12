#!/bin/bash
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  ./runqueries-$file-bfs-sorted.sh > ../../../outputs/ranked/louds/nonFixedQueue/$file.txt
done