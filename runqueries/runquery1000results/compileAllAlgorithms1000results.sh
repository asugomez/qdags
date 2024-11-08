#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# original join
echo "Compiling original join"
cd "$script_dir/original/"
./compileTraditionalJoin.sh

# compile partial results
# louds
echo "Compiling partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./compilePartialLoudsBacktracking.sh
echo "Compiling partial louds optimalOrder"
cd "$script_dir/partial/louds/optimalOrder/"
./compilePartialLoudsNonFixedQueue.sh

# dfuds
echo "Compiling partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktracking.sh
echo "Compiling partial dfuds optimalOrder"
cd "$script_dir/partial/dfuds/optimalOrder/"
./compilePartialDfudsNonFixedQueue.sh

# compile ranked results
# louds
echo "Compiling ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktracking.sh
echo "Compiling ranked louds optimalOrder"
cd "$script_dir/ranked/louds/optimalOrder/"
./compileRankedLoudsNonFixedQueue.sh
# dfuds
echo "Compiling ranked dfuds optimalOrder"
cd "$script_dir/ranked/dfuds/optimalOrder/"
./compileRankedDfudsNonFixedQueue.sh
echo "Compiling ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsBacktracking.sh


## lazy qdags
echo "Compiling lazy join"
cd "$script_dir/lqdags/"
./compileLazyJoin.sh
