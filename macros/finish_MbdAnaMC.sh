#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [generator]"
    exit
fi
source setup_env.sh

generator=$1

dir="${MBD_CENTRALITY_CALIB_PATH}/output/${generator}_magoff/mbdana/"

echo $dir

newfile="${dir}mbd_hist_${generator}_2025.root"

readfiles="${dir}mbd_ana_hist*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 

rm ${readfiles}

newfile="${dir}mbd_tree_${generator}_2025.root"

readfiles="${dir}mbd_ana_tree*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 

rm ${readfiles}

root -q -l "calibration_ana/get_n_events.C(\"${newfile}\")"
