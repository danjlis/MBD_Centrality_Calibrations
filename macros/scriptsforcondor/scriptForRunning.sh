#!/usr/bin/bash                                                                 
runnumber=$1

export CENTRALITY_GAINCALIB=/sphenix/user/dlis/Projects/install/share/centrality/gainfile.calib
export CENTRALITY_BBC_TQ_T0CALIB=/sphenix/user/dlis/Projects/install/share/centrality/bbc_tq_t0.calib

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/Fun4All_CentralityReco.C(${runnumber})"

echo "JOB COMPLETE!"
