#!/bin/bash

###new version###
#argument: <lumi> <systematic> <minseed> <maxseed>

lumi=$1
syst=$2

for i in $(seq $3 $4)
do
    pathstub=`expr $i % 1000`
    logpath=/afs/cern.ch/work/j/joshmt/private/logs/NoEdge/${pathstub}
    if [ ! -d "$logpath" ]; then
	mkdir $logpath
    fi
    bsub -o ${logpath}/EdgeFit_batch_unbinned_${lumi}_${i}_${syst}.log -q 1nh EdgeFit_batch_unbinned_noedge.sh $lumi $syst $i
    #echo "bsub -o ${logpath}/EdgeFit_batch_unbinned_${lumi}_${i}_${syst}.log -q 1nh EdgeFit_batch_unbinned_noedge.sh $lumi $syst $i"

done


