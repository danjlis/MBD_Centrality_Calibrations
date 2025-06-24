#!/bin/bash
if (( $# < 1 )); then
    submitfile=submit_file.txt
else
    submitfile=$1
fi



readarray -t a < $submitfile
njobs=$(cat $submitfile | wc -l)

source setup_env.sh

retDir=${PWD}

cFile=draw_dlis.job
DATE=`date +%Y%m%d`
TIME=`date +%H%M%S`
cDir=condorDir/condor_"$DATE"_"$TIME"

mkdir -p $cDir

cFile2=${cFile%.job}
cFile2="$cFile2"_"$DATE"_"$TIME".job
cp ./job_files/$cFile $cFile2

sed -i -e "s@INITDIR@$PWD/$cDir@g" $cFile2
sed -i -e "s@RUNNUMBER@$runnumber@g" $cFile2

if (( njobs == 1 )); then
    sed -i -e "s@RUNNUMBER@$runnumber@g" $cFile2
else
    sed -i -e "s@RUNNUMBER@0@g" $cFile2
fi
echo "Queue ${njobs}" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForDraw.sh $cDir
cp setup_env.sh $cDir
cp $submitfile $cDir/submit_file.txt

cd $cDir
condor_submit $cFile2
cd $retDir
