#!/bin/bash
for file in p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do # TODO: j3, j4, p2, p3
  ./createExpPri.sh ../../runqueries/original/runqueries-$file-bfs-sorted.sh
done