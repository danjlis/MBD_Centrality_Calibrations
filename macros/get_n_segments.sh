#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage get_n_segments.sh [runnumber]"
    exit
fi

runnumber=${1}

aligned=$(ls /sphenix/lustre01/sphnxpro/commissioning/aligned/*${runnumber}* 2> /dev/null | wc -l)

aligned_prdf=$(ls /sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/*${runnumber}* 2> /dev/null | wc -l)


if (( aligned > aligned_prdf ))
then
    echo $aligned
else
    echo $aligned_prdf
fi

