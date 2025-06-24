#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage submit.sh [runnumber]"
    exit
fi


source setup_env.sh

listfile=$1
runnumber=0
phase=0
if [ ! -f $listfile ]; then
    echo "Nothing"
    exit 0;
fi

mkdir -p ${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/mbdana/
mkdir -p ${MBD_CENTRALITY_CALIB_PATH}/output/run${runnumber}/plots/

cFile=mbdana24_dlis.job
DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=${cFile%.job}
cFile2="$cFile2"_"$DATE"_"$TIME".job

cp ./job_files/$cFile $cFile2

sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2
sed -i -e "s@RUN@$runnumber@g" $cFile2
sed -i -e "s@PHASE@$phase@g" $cFile2

echo "Queue segment from ${listfile}" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForMbdAna24.sh $cDir
cp setup_env.sh $cDir

cd $cDir
condor_submit $cFile2
cd $retDir
