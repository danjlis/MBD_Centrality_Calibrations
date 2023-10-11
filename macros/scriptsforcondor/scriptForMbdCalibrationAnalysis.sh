#!/usr/bin/bash                                                                 
id=$1
runnumber=RUN

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dlis/Projects/install/

source setup_env.sh

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/fun4all/Fun4All_MbdCalibrationAnalysis.C(${runnumber},${id})"

echo "JOB COMPLETE!"
