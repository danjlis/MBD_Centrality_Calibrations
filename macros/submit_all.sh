#!/bin/bash

submitfile=submit_file.txt

if (( $# < 2 )); then
    exit
else
    submitfile=$1
fi

readarray -t a < $submitfile

total_jobs_needed=0
for i in "${a[@]}"
do
    echo "Run ${i} ..."
    bash submit_MbdAna24.sh ${i} 0 ${2}
    sleep 2
    echo "Done."
    
done

echo "Is this bad?"


