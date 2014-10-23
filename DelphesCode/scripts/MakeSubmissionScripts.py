#given a path on eos like
### /eos/cms/store/group/phys_higgs/upgrade/PhaseI/Configuration0/NoPileUp/
#write a .sh file for each subdirectory
#the .sh file should have:
# XXXOUTPATHXXX set to the 'PhaseI/Configuration0/NoPileup/' part of the input path (discarding the rest)
# XXXSAMPLEIDXXX should be set to the name of the subdirectory of the input path
# the sh file should be named XXXSAMPLEIDXXX.sh
# XXXEOSPATHXXX is exactly the input path

#the input .sh script is SubmissionTemplate.sh

def howManyJobs(sample):
  if 'Bj-4p-600' in sample:
    return 6
  elif 'Bj-4p-1100' in sample:
    return 3
  elif 'Bj-4p-300' in sample:
    return 12
  elif 'Bj-4p-0-' in sample:
    return 3
  elif 'Bj' in sample:
    return 2
  elif 'B-' in sample:
    return 2
  elif 'tj-4p-0-600' in sample:
    return 6
  elif 'tj-4p-500' in sample:
    return 4
  elif 'tj-4p-2500' in sample:
    return 3
  elif 'tj' in sample:
    return 2
  elif 'tt-4p-0-' in sample:
    return 4
  elif 'tt-' in sample:
    return 3
  elif 'tB-4p-0-5' in sample:
    return 2
  elif 'tB-4p-1500' in sample:
    return 2
  elif 'tB-4p-2200' in sample:
    return 3
  elif 'tB-4p-500' in sample:
    return 2
  elif 'tB-4p-900' in sample:
    return 4
  elif 'LLB-' in sample:
    return 2
  elif 'LL-4p-100-200-' in sample:
    return 5
  elif 'LL-' in sample:
    return 2

  return 1

import sys,shutil,os,commands
if len(sys.argv) >= 2:
  eospath = sys.argv[1]
else:
  print 'problem -- did not get a path as input'
  sys.exit(1)

script_outdir='submissionScripts/'
if not os.path.exists(script_outdir):
  os.makedirs(script_outdir)

infile=file('.MakeSubmissionScript.tmp','r')
inlist=infile.read().splitlines()
for subdir in inlist:
  #subdir is the sampleid
  #need to decide how many jobs to split into
  njobs = howManyJobs(subdir)
  for ijob in range(1,njobs+1):
    #make subdir.sh file
    scriptfile = script_outdir+subdir+str(ijob)+'.sh'
    print scriptfile
    shutil.copyfile('SubmissionTemplate.sh',scriptfile)
    #we need to find the relevant part of the eospath and put it in outpath
    outpath = eospath[eospath.index('upgrade/')+8:] #yes, this is rather fragile
    cmd1 = 'perl -e "s|XXXOUTPATHXXX|%s|g;" -pi %s' % (outpath,scriptfile)
    (status,out)=commands.getstatusoutput(cmd1) #should, like, test the status or something
    cmd2 = 'perl -e "s|XXXSAMPLEIDXXX|%s|g;" -pi %s' % (subdir,scriptfile)
    (status,out)=commands.getstatusoutput(cmd2) #should, like, test the status or something
    cmd3 = 'perl -e "s|XXXEOSPATHXXX|%s|g;" -pi %s' % (eospath,scriptfile)
    (status,out)=commands.getstatusoutput(cmd3) #should, like, test the status or something
    cmd4 = 'perl -e "s|IJOB|%d|g;" -pi %s' % (ijob,scriptfile)
    (status,out)=commands.getstatusoutput(cmd4) #should, like, test the status or something
    cmd5 = 'perl -e "s|NJOB|%d|g;" -pi %s' % (njobs,scriptfile)
    (status,out)=commands.getstatusoutput(cmd5) #should, like, test the status or something

