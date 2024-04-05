#!/bin/bash

# acces to this file:
# chmod +x ./runqueries/partial/louds/fixedQueue/runqueries-j3-bfs-sorted.sh
# Compilar j3.cpp
#//-std=c++11 -O3 -O0 -DNDEBUG -I /opt/homebrew/include -L /opt/homebrew/lib -lsdsl -ldivsufsort -ldivsufsort64
g++ -std=c++11 -O3 -O0 -DNDEBUG -I /opt/homebrew/include -L /opt/homebrew/lib -lsdsl -ldivsufsort -ldivsufsort64 ../queries/test.cpp -o test

if [ $? -eq 0 ]; then
  ./test ../data/miniTests/grid32R ../data/miniTests/grid32Syzt2 ../data/miniTests/qdagT > ./outputs/test.txt
else
    echo "Error durante la compilación."
fi