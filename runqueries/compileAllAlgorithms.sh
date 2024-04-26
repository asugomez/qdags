#!/bin/bash
script_dir=$(dirname "$BASH_SOURCE")
echo script_dire is $script_dir
# traditional join
cd "$script_dir/all/"
./compileTraditionalJoin.sh

# compile partial results
# louds
cd "$script_dir/partial/louds/backtracking/"
./compilePartialLoudsBacktracking.sh
cd "$script_dir/partial/louds/nonFixedQueue/"
./compilePartialLoudsNonFixedQueue.sh
# dfuds
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktracking.sh
cd "$script_dir/partial/dfuds/nonFixedQueue/"
./compilePartialDfudsNonFixedQueue.sh

# compile ranked results
# louds
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktracking.sh
cd "$script_dir/ranked/louds/nonFixedQueue/"
./compileRankedLoudsNonFixedQueue.sh
# dfuds
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./compileRankedDfudsBacktracking.sh
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsNonFixedQueue.sh
