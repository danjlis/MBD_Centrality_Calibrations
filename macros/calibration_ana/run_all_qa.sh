#!/bin/bash

readarray -t a < runlist.txt

for i in "${a[@]}"
do
    if (( i == 21813 )); then
	continue;
    fi
    echo -n "Run ${i} ..."
    root -l -q "QA_centrality.C(${i})" >> log3.txt
    echo "Done"
done

echo "Done"


