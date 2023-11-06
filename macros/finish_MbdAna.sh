#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [runnumber]"
    exit
fi
source setup_env.sh

runnumber=$1

dir="${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/mbdana/"

echo $dir

newfile="${dir}mbd_trees_${runnumber}.root"

readfiles="${dir}mbd_ana_tree*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 

rm ${readfiles}

root -q -l "calibration_ana/get_n_events.C(\"${newfile}\")"
