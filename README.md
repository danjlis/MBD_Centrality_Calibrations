# FIRST!!#

$ bash setup_dirs.sh

then edit ./macros/setup_env.sh to be the base of this repo


# A couple thigns you need to build before hand #



make a build directory and go there:

$ mkdir build

in build directory

$ mkdir CaloWaveformSimulatorv2
$ cd CaloWaveformSimulatorv2

$../../CaloTriggerEmulator/CaloWaveformSimulatorv2/autogen.sh --prefix=$MYINSTALL
$ make install

in build directory

$ mkdir CaloTriggerEmulator
$ cd CaloTriggerEmulator

$ ../../CaloTriggerEmulator/CaloTriggerEmulator/autogen.sh --prefix=$MYINSTALL
$ make install

in build directory

$ mkdir CaloTriggerAnalysis
$ cd CaloTriggerAnalysis

$ ../../CaloTriggerEmulator/CaloTriggerAnalysis/autogen.sh --prefix=$MYINSTALL
$ make install

in build directory

$ mkdir qacentrality
$ cd qacentrality

$ ../../src/autogen.sh --prefix=$MYINSTALL
$ make install

# This will make everything installed all right (if Centrality Module is in the main coresoftware repo)


MBD Calibrations


$ submit_waveform.sh [runnumber][# of segments]
$ submit_MbdCalibrationAnalysis.sh [runnumber][# of segments]

After these are done

$ submit_singalProcessing.sh [runnumber][# of segments]

$ finish_signals.sh [runnumber]
$ finish_MbdCalibrationAnalysis.sh [runnumber]


Centrality Calibrations

$ submit_centrality.sh [runnumber][# of segments]
$ finish_centrality.sh [runnumber]

These will make your big root files in output/run[runnubmer]