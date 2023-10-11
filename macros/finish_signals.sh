#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage finish.sh [runnumber]"
    exit
fi

runnumber=$1

dir="/sphenix/user/dlis/Projects/signal_processing/output/run${runnumber}/"

echo $dir

newfile="${dir}processedsignal_${runnumber}.root"

readfiles="${dir}signal_*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"

hadd -f ${newfile} ${readfiles} 
