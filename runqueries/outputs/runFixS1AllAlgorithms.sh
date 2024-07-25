#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# partial
## louds
echo "Fixing partial louds backtracking"
./fixS1.sh "$script_dir/partial/louds/backtracking/"
echo "Fixing partial louds nonFixedQueue"
./fixS1.sh "$script_dir/partial/louds/nonFixedQueue/"

## dfuds
echo "Fixing partial dfuds backtracking"
./fixS1.sh "$script_dir/partial/dfuds/backtracking/"
echo "Fixing partial dfuds nonFixedQueue"
./fixS1.sh "$script_dir/partial/dfuds/nonFixedQueue/"

# ranked
## louds
echo "Fixing ranked louds backtracking"
./fixS1.sh "$script_dir/ranked/louds/backtracking/"
echo "Fixing ranked louds nonFixedQueue"
./fixS1.sh "$script_dir/ranked/louds/nonFixedQueue/"

## dfuds
echo "Fixing ranked dfuds nonFixedQueue"
./fixS1.sh "$script_dir/ranked/dfuds/nonFixedQueue/"
echo "Fixing ranked dfuds backtracking"
./fixS1.sh "$script_dir/ranked/dfuds/backtracking/"

## traditional
echo "Fixing traditional join"
./fixS1.sh "$script_dir/all/"

## lazy qdags
echo "Fixing lazy join"
./fixS1.sh "$script_dir/lqdags/"