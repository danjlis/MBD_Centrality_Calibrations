#!/bin/bash

readarray -t a < submit_file.txt

total_jobs_needed=0
for i in "${a[@]}"
do

    segments=$(get_n_segments.sh ${i})

    if [ "${segments}" -gt "5" ]
    then    
	echo ${i} >> good_runs_with_mbd.txt
	echo "Run ${i} ..."
	bash submit_MbdAna.sh ${i}
	sleep 2
	echo "Done."
    fi
done

echo "Is this bad?"


