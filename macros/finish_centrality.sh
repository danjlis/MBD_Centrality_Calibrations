#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [runnumber]"
    exit
fi

source setup_env.sh

runnumber=$1

dir="${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/centrality/"

echo $dir

newfile="${dir}trees_${runnumber}.root"

readfiles="${dir}centrality_reco_tree_*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 
