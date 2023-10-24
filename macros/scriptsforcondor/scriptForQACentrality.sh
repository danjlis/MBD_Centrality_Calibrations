#!/usr/bin/bash                                                                 
id=$1
runnumber=$2

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dlis/Projects/install/

source setup_env.sh

readarray -t a < submit_file.txt

if (( runnumber == 0 )); then
    runnumber=${a[id]} 
fi

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/calibration_ana/QA_centrality.C(${runnumber})"

echo "JOB COMPLETE!"
