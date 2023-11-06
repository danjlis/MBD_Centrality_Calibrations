#!/bin/bash

readarray -t a < good_runs_with_mbd.txt

for i in "${a[@]}"
do

    bash finish_MbdAna.sh ${i}
    
done

echo "done."


