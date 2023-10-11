#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [runnumber]"
    exit
fi

runnumber=$1

dir="/sphenix/user/dlis/Projects/centrality/output/run${runnumber}/"

echo $dir

newfile="${dir}trees_${runnumber}.root"

readfiles="${dir}centrality_reco_tree_*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 
