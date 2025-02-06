#!/bin/bash
for file in j4; do # TODO: j3, j4, p2, p3
  ./createExpPri-j4.sh ../../runqueries/original/runqueries-$file-bfs-sorted.sh
done
