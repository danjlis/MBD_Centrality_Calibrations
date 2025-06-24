#!/usr/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

hostname

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`

if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
    cd $_CONDOR_SCRATCH_DIR
    rsync -a $this_dir/* .
else
   echo condor scratch NOT set
   exit -1
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/dlis/Projects/install/
export SPHENIX=$MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
export LD_LIBRARY_PATH=$MYINSTALL:$LD_LIBRARY_PATH


id=$1
refrun=$2
runnumber=0

source setup_env.sh

readarray -t a < submit_file.txt

runnumber=${a[id]} 

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/calibration_ana/mainreferenceQA.C(${runnumber})"

echo "JOB COMPLETE!"
