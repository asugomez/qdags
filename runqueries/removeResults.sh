#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# partial
## louds
echo "rm partial louds backtracking"
cd "$script_dir/outputs/partial/louds/backtracking/"
rm *
echo "rm partial louds nonFixedQueue"
cd "$script_dir/outputs/partial/louds/nonFixedQueue/"
rm *

## dfuds
echo "rm partial dfuds backtracking"
cd "$script_dir/outputs/partial/dfuds/backtracking/"
rm *
echo "rm partial dfuds nonFixedQueue"
cd "$script_dir/outputs/partial/dfuds/nonFixedQueue/"
rm *

# ranked
## louds
echo "rm ranked louds backtracking"
cd "$script_dir/outputs/ranked/louds/backtracking/"
rm *
echo "rm ranked louds nonFixedQueue"
cd "$script_dir/outputs/ranked/louds/nonFixedQueue/"
rm *

## dfuds
echo "rm ranked dfuds nonFixedQueue"
cd "$script_dir/outputs/ranked/dfuds/nonFixedQueue/"
rm *
echo "rm ranked dfuds backtracking"
cd "$script_dir/outputs/ranked/dfuds/backtracking/"
rm *

## traditional
echo "rm traditional join"
cd "$script_dir/outputs/all/"
rm *


## lazy qdags
echo "rm lazy qdags join"
cd "$script_dir/outputs/lqdags/"
rm *
