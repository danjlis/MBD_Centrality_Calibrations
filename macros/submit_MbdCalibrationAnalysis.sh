#!/bin/bash

if (( $# < 2 )); then
    echo "Illegal number of parameters"
    echo "Usage submit.sh [runnumber] [# of segments]"
    exit
fi

source setup_env.sh

runnumber=$1
retDir=${PWD}

mkdir -p  ${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/mbdcalibana/
mkdir -p  ${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/plots/
cFile=mbdcalibana_dlis.job
DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=${cFile%.job}
cFile2="$cFile2"_"$DATE"_"$TIME".job
cp ./job_files/$cFile $cFile2


sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2

echo "Queue $2" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForMbdCalibrationAnalysis.sh $cDir
cp setup_env.sh $cDir
sed -i -e "s@RUN@$runnumber@g" $cDir/scriptForMbdCalibrationAnalysis.sh

cd $cDir
condor_submit $cFile2
cd $retDir
