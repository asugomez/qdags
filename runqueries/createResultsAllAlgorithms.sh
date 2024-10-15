#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"


### original
#echo "Results original join"
#cd "$script_dir/original/"
#./createResultsSummary.sh

# partial
## louds
echo "Results partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./createResultsSummary.sh
echo "Results partial louds optimalOrder"
cd "$script_dir/partial/louds/optimalOrder/"
./createResultsSummary.sh

## dfuds
echo "Results partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./createResultsSummary.sh
#echo "Results partial dfuds optimalOrder"
#cd "$script_dir/partial/dfuds/optimalOrder/"
#./createResultsSummary.sh
#
## ranked
### louds
#echo "Results ranked louds backtracking"
#cd "$script_dir/ranked/louds/backtracking/"
#./createResultsSummary.sh
#echo "Results ranked louds optimalOrder"
#cd "$script_dir/ranked/louds/optimalOrder/"
#./createResultsSummary.sh
#
### dfuds
#echo "Results ranked dfuds optimalOrder"
#cd "$script_dir/ranked/dfuds/optimalOrder/"
#./createResultsSummary.sh
#echo "Results ranked dfuds backtracking"
#cd "$script_dir/ranked/dfuds/backtracking/"
#./createResultsSummary.sh

#
### lazy qdags
#echo "Results lazy join"
#cd "$script_dir/lqdags/"
#./createResultsSummary.sh