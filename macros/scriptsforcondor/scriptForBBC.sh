#!/usr/bin/bash                                                                 
runnumber=$1

export BBCCALIB=/sphenix/user/dlis/Projects/centrality/calib/

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/Fun4All_CentralityRecoMBD.C(${runnumber})"

echo "JOB COMPLETE!"
