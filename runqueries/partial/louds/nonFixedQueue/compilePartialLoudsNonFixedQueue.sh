#!/bin/bash
# acces to this file:
# chmod +x ./runqueries/partial/louds/backtracking/runqueries-j3-bfs-sorted.sh
# Compilar j3.cpp
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  g++ -std=c++11 -O3 -O0 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 "../../../../queries/partial/louds/nonFixedQueue/$file.cpp" -o $file
done
