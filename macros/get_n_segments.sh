#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage get_n_segments.sh [runnumber]"
    exit
fi

aligned=$(ls /sphenix/lustre01/sphnxpro/commissioning/aligned/*${1}* | wc -l)

aligned_prdf=$(ls /sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/*${1}* | wc -l)

if (( aligned > 0 ))
then
    echo "Run ${1} has ${aligned} segments in sphenix/lustre01/sphnxpro/commissioning/aligned/"
fi

if (( aligned_prdf > 0 ))
then
    echo "Run ${1} has ${aligned_prdf} segments in sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/"
fi

if (( aligned_prdf == 0 && aligned == 0))
then
    echo "Run ${1} has NOTHING"
fi



