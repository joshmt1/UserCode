#!/bin/bash

#hack for inserting lines for lhapdffull libraries into CMSSW.sh
#usage:
#crab -create
#./crabHack.sh crabdir
#crab -submit

if [ $# -ne 1 ]
then
  echo "Usage: ./crabHack.sh crabdir"
  exit 1
fi

#take crab dir as argument
echo $1

#search for SCRAMRT_LSB_JOBNAME
#the lines need to be inserted *after* the first 'fi' after that string
thelinenum=`awk '/SCRAMRT_LSB_JOBNAME/ {print NR; exit;}' $1/job/CMSSW.sh`

#must use double quotes for the dollar sign to work!
linenum2=`awk "/fi/ {if (NR> $thelinenum) {print NR; exit;}}" $1/job/CMSSW.sh`

linenum2=`expr $linenum2 + 1`

#insert:
#scramv1 setup lhapdffull
#scramv1 b
sed "${linenum2}i\scramv1 setup lhapdffull\nscramv1 b" $1/job/CMSSW.sh > $1/job/CMSSW.new.sh

#just to show the user the results
grep -B7 -A2 lhapdffull $1/job/CMSSW.new.sh

mv $1/job/CMSSW.new.sh $1/job/CMSSW.sh

echo " ==================== "
echo "Ready for crab -submit"
