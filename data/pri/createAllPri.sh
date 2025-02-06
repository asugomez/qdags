#!/bin/bash
for file in p2; do # TODO: j3, j4, p2, p3
  ./createExpPri.sh ../../runqueries/original/runqueries-$file-bfs-sorted.sh
done
