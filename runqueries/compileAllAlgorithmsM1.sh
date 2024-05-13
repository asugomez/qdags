#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# compile partial results
# louds
cd "$script_dir/partial/louds/backtracking/"
./compilePartialLoudsBacktrackingM1.sh
cd "$script_dir/partial/louds/nonFixedQueue/"
./compilePartialLoudsNonFixedQueueM1.sh
# dfuds
cd "$script_dir/partial/dfuds/backtracking/"
./compilePartialDfudsBacktrackingM1.sh
cd "$script_dir/partial/dfuds/nonFixedQueue/"
./compilePartialDfudsNonFixedQueueM1.sh

# compile ranked results
# louds
cd "$script_dir/ranked/louds/backtracking/"
./compileRankedLoudsBacktrackingM1.sh
cd "$script_dir/ranked/louds/nonFixedQueue/"
./compileRankedLoudsNonFixedQueueM1.sh
# dfuds
cd "$script_dir/ranked/dfuds/nonFixedQueue/"
./compileRankedDfudsBacktrackingM1.sh
cd "$script_dir/ranked/dfuds/backtracking/"
./compileRankedDfudsNonFixedQueueM1.sh

# traditional join
cd "$script_dir/all/"
./compileTraditionalJoinM1.sh
