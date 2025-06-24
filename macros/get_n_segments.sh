#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage get_n_segments.sh [runnumber]"
    exit
fi

runnumber=${1}

mine=$(ls /sphenix/user/dlis/Projects/zdc_fix/output/*${runnumber}* 2> /dev/null | wc -l)

#aligned=$(ls /sphenix/lustre01/sphnxpro/commissioning/mbd/beam/beam_seb18-000${runnumber}* 2> /dev/null | wc -l)

aligned_prdf=$(ls /sphenix/lustre01/sphnxpro/commissioning/aligned_prdf/*${runnumber}* 2> /dev/null | wc -l)



if (( mine > aligned_prdf ))
then
    echo $mine
else
    echo $aligned_prdf
fi

