#!/bin/bash

submitfile=$1

readarray -t a < $submitfile

for i in "${a[@]}"
do
    
    bash finish_MbdAna.sh ${i}
    
done

echo "done."


