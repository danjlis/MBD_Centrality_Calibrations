#!/usr/bin/bash                                                                 
id=$1
runnumber=$2

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dlis/Projects/install/

source setup_env.sh

readarray -t a < submit_file.txt


runnumber=${a[id]} 

echo ${runnumber}

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/drawing/Draw_QA_Centrality.C(${runnumber})"

echo "JOB COMPLETE!"
