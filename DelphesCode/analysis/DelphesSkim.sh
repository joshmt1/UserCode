#! /bin/bash 

if [ ${#1} = 0 ]; then
echo "Requires a directory as an argument"
exit
fi

#don't skim signal
mylist=`ls -1 $1/reduced*.root | grep -v susy | grep -v natural`

#MLL TIGHTHT
#export SKIMOPTION='TIGHTHT MLL'
#export SKIMOPTION='MLL'
export SKIMOPTION='TIGHTHT'

for j in $mylist;
do 

echo $j
root -b -l -q "DelphesSkim.C+(\"$j\")"

#ps aux | grep joshmt | grep 'root.exe' | grep -v grep | wc -l

done

#ps aux | grep joshmt | grep -v grep | grep -c 'root\.exe'
