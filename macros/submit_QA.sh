#!/bin/bash

if (( $# < 1 )); then
    readarray -t a < good_runs_with_mbd.txt
    njobs=$(cat good_runs_with_mbd.txt | wc -l)
else
    njobs=1
    runnumber=$1
fi

source setup_env.sh

retDir=${PWD}

cFile=qacentrality_dlis.job
DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=${cFile%.job}
cFile2="$cFile2"_"$DATE"_"$TIME".job
cp ./job_files/$cFile $cFile2

sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2

if (( njobs == 1 )); then
    sed -i -e "s@ARGS@$runnumber@g" $cFile2
else
    sed -i -e "s@ARGS@0@g" $cFile2
fi
echo "Queue ${njobs}" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForQACentrality.sh $cDir
cp setup_env.sh $cDir
cp submit_file.txt $cDir

cd $cDir
condor_submit $cFile2
cd $retDir
