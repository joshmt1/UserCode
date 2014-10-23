#!/bin/sh
echo "~~~ Started"
workdir=$PWD
exedir=/afs/cern.ch/work/j/joshmt/private/Upgrade/Delphes
outputdir=/eos/cms/store/user/joshmt/upgrade/v10/XXXOUTPATHXXX
outfilename=reducedTree.def.XXXSAMPLEIDXXX_IJOB_NJOB.root

eoscmd=/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select

source /afs/cern.ch/sw/lcg/external/gcc/4.7.2/x86_64-slc6-gcc47-opt/setup.sh
cd /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc47-opt/root
source bin/thisroot.sh
cd $exedir
echo "FlatTree starting"
./FlatTree "root://eoscms.cern.ch/XXXEOSPATHXXX/XXXSAMPLEIDXXX/*.root" -O ${workdir}/${outfilename} -N IJOB NJOB 
ls -l ${workdir}/*.root
${eoscmd} mkdir -p $outputdir
${eoscmd} cp ${workdir}/${outfilename} ${outputdir}/
