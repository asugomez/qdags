#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

## original
echo "Running original join"
cd "$script_dir/original/"
./runTraditionalJoin.sh


# ranked
## louds
echo "Running ranked louds backtracking"
cd "$script_dir/ranked/louds/backtracking/"
./runRankedLoudsBacktracking.sh
echo "Running ranked louds optimalOrder"
cd "$script_dir/ranked/louds/optimalOrder/"
./runRankedLoudsNonFixed.sh

## dfuds
echo "Running ranked dfuds optimalOrder"
cd "$script_dir/ranked/dfuds/optimalOrder/"
./runRankedDfudsNonFixed.sh
echo "Running ranked dfuds backtracking"
cd "$script_dir/ranked/dfuds/backtracking/"
./runRankedDfudsBacktracking.sh

# partial
# louds
echo "Running partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./runPartialLoudsBacktracking.sh
echo "Running partial louds optimalOrder"
cd "$script_dir/partial/louds/optimalOrder/"
./runPartialLoudsNonFixed.sh

## dfuds
echo "Running partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./runPartialDfudsBacktracking.sh
echo "Running partial dfuds optimalOrder"
cd "$script_dir/partial/dfuds/optimalOrder/"
./runPartialDfudsNonFixed.sh

## lazy qdags
#echo "Running lazy join"
#cd "$script_dir/lqdags/"
#./runLazyJoin.sh
