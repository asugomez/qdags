#!/bin/bash

for file in j4; do #j3 p2 p3 p4 s1 s2 s3 s4 t2 t3 t4 ti2 ti3 ti4 tr1 tr2; do
  archivo="./runqueries-$file-bfs-sorted.sh"
  newName="./runqueries-$file-bfs-sorted-fun-1.sh" # change number to 0,1,2...
  while IFS= read -r linea || [[ -n "$linea" ]]; do
      echo "$linea 1" # change number to 0,1,2...
  done < "$archivo" > "$archivo.tmp"
  mv "$archivo.tmp" "$newName"
done
