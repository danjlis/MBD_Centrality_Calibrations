#!/usr/bin/bash                                                                 
id=$1
runnumber=RUN

source /opt/sphenix/core/bin/sphenix_setup.sh -n
source /opt/sphenix/core/bin/setup_local.sh /sphenix/user/dlis/Projects/install/

source setup_env.sh

export CENTRALITY_GAINCALIB=/sphenix/user/dlis/Projects/install/share/centrality/gaincalib.calib
export CENTRALITY_BBC_TQ_T0CALIB=/sphenix/user/dlis/Projects/install/share/centrality/bbc_tq_t0.calib

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/fuun4all/Fun4All_CentralityRecoQA.C(${runnumber},${id})"

echo "JOB COMPLETE!"
