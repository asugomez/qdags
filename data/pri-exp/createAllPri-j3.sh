#!/bin/bash
for file in j3; do # TODO: j3, j4, p2, p3
  ./createExpPri-j3.sh ../../runqueries/original/runqueries-$file-bfs-sorted.sh
done
