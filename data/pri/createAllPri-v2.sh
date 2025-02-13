#!/bin/bash
for file in p4; do # TODO: j3, j4, p2, p3
  ./createExpPri-v2.sh ../../runqueries/original/runqueries-$file-bfs-sorted.sh
done
