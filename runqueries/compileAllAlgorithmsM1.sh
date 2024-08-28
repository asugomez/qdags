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
echo "Compiling partial louds nonFixedQueue"
cd "$script_dir/partial/louds/nonFixedQueue/"
./compilePartialLoudsNonFixedQueueM1.sh
# dfuds
echo "Compiling partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktrackingM1.sh
echo "Compiling partial dfuds nonFixedQueue"
cd "$script_dir/partial/dfuds/nonFixedQueue/"
./compilePartialDfudsNonFixedQueueM1.sh

# compile ranked results
# louds
echo "Compiling ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktrackingM1.sh
echo "Compiling ranked louds nonFixedQueue"
cd "$script_dir/ranked/louds/nonFixedQueue/"
./compileRankedLoudsNonFixedQueueM1.sh
# dfuds
echo "Compiling ranked dfuds nonFixedQueue"
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./compileRankedDfudsNonFixedQueueM1.sh
echo "Compiling ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsBacktrackingM1.sh


## lazy qdags
echo "Compiling lazy join"
cd "$script_dir/lqdags/"
./compileLazyJoinM1.sh