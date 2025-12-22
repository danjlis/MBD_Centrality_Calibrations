# FIRST!!#

$ bash setup_dirs.sh

then edit ./macros/setup_env.sh to be the base of this repo

# A couple things you need to build before hand #


make a build directory and go there:

$ mkdir build

in build directory

$ mkdir MbdCalibrationAnalysis
$ cd MbdCalibrationAnalysis

$../../src/MbdCalibrationAnalysis/autogen.sh --prefix=$MYINSTALL
$ make install

in build directory

$ mkdir CaloTriggerEmulator
$ cd CaloTriggerEmulator

$ ../../CaloTriggerEmulator/CaloTriggerEmulator/autogen.sh --prefix=$MYINSTALL
$ make install

# This will make everything installed all right (if Centrality Module is in the main coresoftware repo)


These will make your big root files in output/run[runnubmer]
