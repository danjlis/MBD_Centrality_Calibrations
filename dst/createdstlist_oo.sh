#! /usr/bin/env bash
# maybe replace cdb and build with arguments
# get all runs in a list

#make the dst files in the run

rm listfile.list

CreateDstList.pl --dataset run3oo --printruns --tag ana536_2025p010_v001 DST_CALOFITTING >> listfile.list

CreateDstList.pl --dataset run3oo  --tag ana536_2025p010_v001 DST_CALOFITTING --list listfile.list



