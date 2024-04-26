#!/bin/bash
script_dir=$(dirname "$BASH_SOURCE")
for file in j3 j4 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  g++ -std=c++11 -O3 -O0 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 "$pwd/../../../../queries/ranked/louds/backtracking/$file.cpp" -o "$script_dir/$file"
done
