#!/bin/sh
echo "~~~ Started"


lumi=$1
syst=$2
seed=$3

#path="DelphesEdge/NoEdge/0"

cd /afs/cern.ch/work/j/joshmt/private/SusyRa2b/CMSSW_5_3_11/src
eval `scramv1 runtime -sh`
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh
cd NtupleTools

commandstringgen="EdgeFit_genone.C(${lumi},${seed})"
root -b -l -q $commandstringgen
commandstringfit="EdgeFit_runoneunbinned_noedge.C(${lumi},${seed},${syst})"
root -b -l -q $commandstringfit
# >& ${path}/run_one_unbinned_${lumi}_${jobnumber}_${syst}.log

