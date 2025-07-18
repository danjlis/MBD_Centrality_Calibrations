#!/bin/bash

if (( $# < 1 )); then
    echo "Illegal number of parameters"
    echo "Usage submit.sh [runnumber]"
    exit
fi


source setup_env.sh

runnumber=$1

retDir=${PWD}

cFile=mbdana24_dlis.job
DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=${cFile%.job}
cFile2="$cFile2"_"$DATE"_"$TIME".job
cp ./job_files/$cFile $cFile2

sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2

echo "Queue ${segments}" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForMbdAna24.sh $cDir
cp setup_env.sh $cDir
sed -i -e "s@RUN@$runnumber@g" $cDir/scriptForMbdAna24.sh

cd $cDir
condor_submit $cFile2
cd $retDir
