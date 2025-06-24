#!/usr/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

hostname

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/dlis/Projects/install/
export SPHENIX=$MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
export LD_LIBRARY_PATH=$MYINSTALL:$LD_LIBRARY_PATH

#source /cvmfs/sphenix.sdcc.bnl.gov/gcc-12.1.0/opt/sphenix/core/bin/setup_local.sh /sphenix/u/sphnxpro/chp/eventcombine/install

if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
    cd $_CONDOR_SCRATCH_DIR
    rsync -a $this_dir/* .
else
   echo condor scratch NOT set
   exit -1
fi

runnumber=$1
refrun=$2

source setup_env.sh

root -b -q "/sphenix/user/dlis/Projects/centrality/macros/calibration_ana/mainQA.C(${runnumber}, ${refrun})"

echo "JOB COMPLETE!"
