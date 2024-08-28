#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# partial
## louds
echo "Fixing partial louds backtracking"
./fixS1.sh "$script_dir/partial/louds/backtracking/"
echo "Fixing partial louds optimalOrder"
./fixS1.sh "$script_dir/partial/louds/optimalOrder/"

## dfuds
echo "Fixing partial dfuds backtracking"
./fixS1.sh "$script_dir/partial/dfuds/backtracking/"
echo "Fixing partial dfuds optimalOrder"
./fixS1.sh "$script_dir/partial/dfuds/optimalOrder/"

# ranked
## louds
echo "Fixing ranked louds backtracking"
./fixS1.sh "$script_dir/ranked/louds/backtracking/"
echo "Fixing ranked louds optimalOrder"
./fixS1.sh "$script_dir/ranked/louds/optimalOrder/"

## dfuds
echo "Fixing ranked dfuds optimalOrder"
./fixS1.sh "$script_dir/ranked/dfuds/optimalOrder/"
echo "Fixing ranked dfuds backtracking"
./fixS1.sh "$script_dir/ranked/dfuds/backtracking/"

## traditional
echo "Fixing traditional join"
./fixS1.sh "$script_dir/all/"

## lazy qdags
echo "Fixing lazy join"
./fixS1.sh "$script_dir/lqdags/"