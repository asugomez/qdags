#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# original join
echo "Compiling original join"
cd "$script_dir/original/"
./compileTraditionalJoinM1.sh

# compile partial results
# louds
echo "Compiling partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./compilePartialLoudsBacktrackingM1.sh
echo "Compiling partial louds optimalOrder"
cd "$script_dir/partial/louds/optimalOrder/"
./compilePartialLoudsNonFixedQueueM1.sh
# dfuds
echo "Compiling partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktrackingM1.sh
echo "Compiling partial dfuds optimalOrder"
cd "$script_dir/partial/dfuds/optimalOrder/"
./compilePartialDfudsNonFixedQueueM1.sh

# compile ranked results
# louds
echo "Compiling ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktrackingM1.sh
echo "Compiling ranked louds optimalOrder"
cd "$script_dir/ranked/louds/optimalOrder/"
./compileRankedLoudsNonFixedQueueM1.sh
# dfuds
echo "Compiling ranked dfuds optimalOrder"
cd "$script_dir/ranked/dfuds/optimalOrder/"
./compileRankedDfudsNonFixedQueueM1.sh
echo "Compiling ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsBacktrackingM1.sh

## lazy qdags
echo "Compiling lazy join"
cd "$script_dir/lqdags/"
./compileLazyJoinM1.sh