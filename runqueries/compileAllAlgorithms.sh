#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# compile partial results
# louds
echo "Compiling partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./compilePartialLoudsBacktracking.sh
echo "Compiling partial louds nonFixedQueue"
cd "$script_dir/partial/louds/nonFixedQueue/"
./compilePartialLoudsNonFixedQueue.sh

# dfuds
echo "Compiling partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktracking.sh
echo "Compiling partial dfuds nonFixedQueue"
cd "$script_dir/partial/dfuds/nonFixedQueue/"
./compilePartialDfudsNonFixedQueue.sh

# compile ranked results
# louds
echo "Compiling ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktracking.sh
echo "Compiling ranked louds nonFixedQueue"
cd "$script_dir/ranked/louds/nonFixedQueue/"
./compileRankedLoudsNonFixedQueue.sh
# dfuds
echo "Compiling ranked dfuds nonFixedQueue"
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./compileRankedDfudsNonFixedQueue.sh
echo "Compiling ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsBacktracking.sh

# traditional join
echo "Compiling traditional join"
cd "$script_dir/all/"
./compileTraditionalJoin.sh
