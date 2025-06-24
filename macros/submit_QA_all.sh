#!/bin/bash
if (( $# < 2 )); then
    submitfile=submit_file.txt
    referencerun=23696
else
    submitfile=$1
    referencerun=$2
fi

readarray -t a < $submitfile
njobs=$(cat $submitfile | wc -l)

echo $runnumber

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
sed -i -e "s@REFRUN@$referencerun@g" $cFile2
sed -i -e "s@RUNNUMBER@0@g" $cFile2
echo "Queue ${njobs}" >> $cFile2

cp $cFile2 $cDir
rm $cFile2
cp scriptsforcondor/scriptForQACentrality.sh $cDir
cp setup_env.sh $cDir
cp $submitfile $cDir/submit_file.txt

cd $cDir
condor_submit $cFile2
cd $retDir
