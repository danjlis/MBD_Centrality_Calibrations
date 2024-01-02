
export MYINSTALL=/sphenix/user/dlis/Projects/install/
source /opt/sphenix/core/bin/sphenix_setup.sh -n new
export SPHENIX=$MYINSTALL
source /opt/sphenix/core/bin/setup_local.sh $MYINSTALL
export LD_LIBRARY_PATH=$MYINSTALL:$LD_LIBRARY_PATH
#ROOT
#source /opt/sphenix/core/root/bin/thisroot.sh

#PYTHIA
#export PYTHIA8PATH=/opt/sphenix/core/pythia8/bin
#export PATH=$PATH:$PYTHIA8PATH

#FASTJET
#export FASTJETPATH=/opt/sphenix/core/fastjet/bin
#export PATH=$PATH:$FASTJETPATH

#SPHENIX
#DEFAULT
#source /opt/sphenix/core/bin/sphenix_setup.sh -n new
#PLAY
#source /opt/sphenix/core/bin/sphenix_setup.sh -n play.3

#source /opt/sphenix/core/bin/setup_local.sh 
#export TESTINSTALL=/sphenix/u/cmcginn/testinstall
#export LD_LIBRARY_PATH=$TESTINSTALL/lib:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/sphenix/u/cmcginn/fun4allbuild/lib:$LD_LIBRARY_PATH

