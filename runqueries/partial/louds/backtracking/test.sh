#!/bin/bash
g++ -std=c++11 -O3 -O0 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 ./queries/partial/louds/backtracking/j3.cpp -o ./j3test1
g++ -std=c++11 -O3 -O0 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 ./../../../../queries/partial/louds/backtracking/j3.cpp -o ./j3test2