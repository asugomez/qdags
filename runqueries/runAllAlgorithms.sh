#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# partial
## louds
echo "Running partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./runPartialLoudsBacktracking.sh
echo "Running partial louds nonFixedQueue"
cd "$script_dir/partial/louds/nonFixedQueue/"
./runPartialLoudsNonFixed.sh

## dfuds
echo "Running partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./runPartialDfudsBacktracking.sh
echo "Running partial dfuds nonFixedQueue"
cd "$script_dir/partial/dfuds/nonFixedQueue/"
./runPartialDfudsNonFixedQueue.sh

# ranked
## louds
echo "Running ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./runRankedLoudsBacktracking.sh
echo "Running ranked louds nonFixedQueue"
cd "$script_dir/ranked/louds/nonFixedQueue/"
./runRankedLoudsNonFixedQueue.sh

## dfuds
echo "Running ranked dfuds nonFixedQueue"
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./runRankedDfudsNonFixedQueue.sh
echo "Running ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./runRankedDfudsBacktracking.sh

## traditional
echo "Running traditional join"
cd "$script_dir/all/"
./runTraditionalJoin.sh
