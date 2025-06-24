#!/usr/bin/bash

export USER="$(id -u -n)"
export LOGNAME=${USER}
export HOME=/sphenix/u/${USER}

this_script=$BASH_SOURCE
this_script=`readlink -f $this_script`
this_dir=`dirname $this_script`


container_script=container_`basename $this_script`
singularity exec -B /home -B /direct/sphenix+u -B /gpfs02 -B /sphenix/u -B /sphenix/lustre01 -B /sphenix/user  -B /sphenix/sim /cvmfs/sphenix.sdcc.bnl.gov/singularity/rhic_sl7.sif ./$container_script $*
echo wrapper script done
