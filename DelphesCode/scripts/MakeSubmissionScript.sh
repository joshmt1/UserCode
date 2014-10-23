###no we don't want this!
#rm -rf submissionScripts/
eoscms ls $1 > .MakeSubmissionScript.tmp
python MakeSubmissionScripts.py $1
rm .MakeSubmissionScript.tmp
chmod +x submissionScripts/*.sh
