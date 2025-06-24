#!/usr/bin/bash                                                                 

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`
hostname

if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]
then
    cd $_CONDOR_SCRATCH_DIR
    rsync -a $this_dir/* .
else
   echo condor scratch NOT set
   exit -1
fi

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.489
export MYINSTALL=/sphenix/user/dlis/Projects/install/
export SPHENIX=$MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
export LD_LIBRARY_PATH=$MYINSTALL:$LD_LIBRARY_PATH

id=$1
phase=$2
source setup_env.sh

root -b -q "Fun4All_MbdAnaDST24.C(\"${id}\", 0)"

cp *out /sphenix/user/dlis/Projects/centrality/macros/logs/
cp *err /sphenix/user/dlis/Projects/centrality/macros/logs/
cp *root /sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/centrality/

echo "JOB COMPLETE!"

