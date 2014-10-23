#!/bin/bash

###new version###
#argument: <input path> <systematics arg>

#example:
#./EdgeFit_batch_unbinned_submit.sh /afs/cern.ch/work/j/joshmt/private/SusyRa2b/CMSSW_5_3_11/src/NtupleTools/DelphesEdge/3000 -1
#2nd argument is the constraint to use in the fit (-1 = no constraint)


inputdir=$1
inputfiles=`ls $inputdir/datasets_*.root`

pattern="_([0-9]+)_([0-9]+).root"

for file in $inputfiles
do
    echo $file
    #bash regex
    [[ $file =~ $pattern ]]
    lumi="${BASH_REMATCH[1]}"
    index="${BASH_REMATCH[2]}"
    #echo $lumi
    #echo $index
    #echo $(dirname $file)
    bsub -o /afs/cern.ch/work/j/joshmt/private/logs/EdgeFit_batch_unbinned_${lumi}_${index}_$2.log -q 1nh EdgeFit_batch_unbinned.sh $lumi $2 $file
    #sleep 1
done
