#!/bin/bash
for file in p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  ./fixNumArgs.sh ./runqueries-$file-bfs-sorted.sh
done