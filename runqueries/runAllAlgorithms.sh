#!/bin/bash
script_dir=$(dirname "$BASH_SOURCE")

# partial
## louds
cd "$script_dir/partial/louds/backtracking/"
./runPartialLoudsBacktracking.sh

cd "$script_dir/partial/louds/nonFixedQueue/"
./runPartialLoudsNonFixed.sh

## dfuds
cd "$script_dir/partial/dfuds/backtracking/"
./runPartialDfudsBacktracking.sh

cd "$script_dir/partial/dfuds/nonFixedQueue/"
./runPartialDfudsNonFixedQueue.sh

# ranked
## louds
cd "$script_dir/ranked/louds/backtracking/"
./runRankedLoudsBacktracking.sh

cd "$script_dir/ranked/louds/nonFixedQueue/"
./runRankedLoudsNonFixedQueue.sh

## dfuds
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./runRankedDfudsNonFixedQueue.sh

cd "$script_dir/ranked/dfuds/backtracking/"
./runRankedDfudsBacktracking.sh
