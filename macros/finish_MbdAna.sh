#!/bin/bash


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

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export MYINSTALL=/sphenix/user/dlis/Projects/install/
export SPHENIX=$MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
export LD_LIBRARY_PATH=$MYINSTALL:$LD_LIBRARY_PATH

source setup_env.sh

runnumber=$1

dir="/sphenix/tg/tg01/commissioning/CaloCalibWG/dlis/centrality/"
hdir="${dir}runs/"

mkdir -p $hdir

echo $dir

newfile="${hdir}/mbd_trees_${runnumber}.root"

readfiles="${dir}/mbd_ana_tree_000${runnumber}*.root"

echo "From: ${readfiles}"
echo "To: ${newfile}"
if [ -f $newfile ]; then
    echo "File exists"
#    exit 0;
fi

hadd -f ${newfile} ${readfiles}

