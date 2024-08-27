#!/bin/bash
# choose one dataset from ./data/original

script_dir=$(dirname "$BASH_SOURCE")
datasets=("$script_dir/original"/*)
numDatasets=${#datasets[@]}

# random number between 0 and numDatasets
index=$((RANDOM % numDatasets))

dataset_selected="${datasets[$index]}"
name_file=$(basename "$dataset_selected")

echo $name_file