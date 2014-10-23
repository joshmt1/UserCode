#!/bin/sh
echo "~~~ Started"

#lumi=3000
#jobnumber=109000
#syst=0.1

lumi=$1
path=$3
syst=$2

#path="DelphesEdge/NoEdge/0"

cd /afs/cern.ch/work/j/joshmt/private/SusyRa2b/CMSSW_5_3_11/src
eval `scramv1 runtime -sh`
source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.18/x86_64-slc5-gcc46-opt/root/bin/thisroot.sh
cd NtupleTools

commandstring="EdgeFit_runoneunbinned.C(${lumi},\"${path}\",${syst})"
root -b -l -q $commandstring
# >& ${path}/run_one_unbinned_${lumi}_${jobnumber}_${syst}.log

