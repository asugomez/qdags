#!/bin/bash
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# partial
## louds
echo "Results partial louds backtracking"
cd "$script_dir/partial/louds/backtracking/"
./createResultsSummary.sh
echo "Results partial louds nonFixedQueue"
cd "$script_dir/partial/louds/nonFixedQueue/"
./createResultsSummary.sh

## dfuds
echo "Results partial dfuds backtracking"
cd "$script_dir/partial/dfuds/backtracking/"
./createResultsSummary.sh
#echo "Results partial dfuds nonFixedQueue"
#cd "$script_dir/partial/dfuds/nonFixedQueue/"
#./createResultsSummary.sh
#
## ranked
### louds
#echo "Results ranked louds backtracking"
#cd "$script_dir/ranked/louds/backtracking/"
#./createResultsSummary.sh
#echo "Results ranked louds nonFixedQueue"
#cd "$script_dir/ranked/louds/nonFixedQueue/"
#./createResultsSummary.sh
#
### dfuds
#echo "Results ranked dfuds nonFixedQueue"
#cd "$script_dir/ranked/dfuds/nonFixedQueue/"
#./createResultsSummary.sh
#echo "Results ranked dfuds backtracking"
#cd "$script_dir/ranked/dfuds/backtracking/"
#./createResultsSummary.sh

### traditional
#echo "Results traditional join"
#cd "$script_dir/all/"
#./createResultsSummary.sh
#
### lazy qdags
#echo "Results lazy join"
#cd "$script_dir/lqdags/"
#./createResultsSummary.sh