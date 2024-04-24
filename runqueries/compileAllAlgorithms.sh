#!/bin/bash
# traditional join
./all/compileTraditionalJoin.sh

# compile partial results
# dfuds
./partial/dfuds/backtracking/compilePartialDfudsBacktracking.sh
./partial/dfuds/nonFixedQueue/compilePartialDfudsNonFixedQueue.sh
# louds
./partial/louds/backtracking/compilePartialLoudsBacktracking.sh
./partial/louds/nonFixedQueue/compilePartialLoudsNonFixedQueue.sh

# compile ranked results
# dfuds
./ranked/dfuds/backtracking/compileRankedDfudsBacktracking.sh
./ranked/dfuds/nonFixedQueue/compileRankedDfudsNonFixedQueue.sh
# louds
./ranked/louds/backtracking/compileRankedLoudsBacktracking.sh
./ranked/louds/nonFixedQueue/compileRankedLoudsNonFixedQueue.sh
