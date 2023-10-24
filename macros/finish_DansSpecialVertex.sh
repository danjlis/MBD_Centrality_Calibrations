Srm#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [runnumber]"
    exit
fi
source setup_env.sh

runnumber=$1

dir="${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/plots/"

echo $dir

newfile="${dir}dans_h_vtx_${runnumber}.root"

readfiles="${dir}dans_vtx_*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 

rm ${readfiles}

