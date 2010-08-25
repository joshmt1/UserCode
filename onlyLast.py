import os, sys, re

###############################################################

#Usage:
# python onlyLast.py <path to input directory on CASTOR>

#now updated to handle output from the latest CRAB with
#file names like this:
#                      myfilename_1_1_aBC.root

###############################################################

#handle input arguments
#this includes the name of the script
##for arg in sys.argv:
##    print arg
    
#sys.argv[1] will have the first argument
tmpfile = '.onlyLast_py_'
mypid = os.getpid()
tmpfile += `mypid`

inputdir = sys.argv[1]
if inputdir[len(inputdir)-1] != '/':
    inputdir+='/'

outputdir = inputdir
outputdir += 'EXTRAS'

mkdircommand = "nsmkdir "
mkdircommand += outputdir

#get a file listing into tmpfile
cmd = 'nsls -l '
cmd += inputdir
cmd += ' | awk \'// {print $5,$9;}\' > '
cmd += tmpfile
#print cmd
os.system(cmd)

#keep track of whether we made the directory or not
alreadymadedir = 0

f = open(tmpfile,'r')

#to deal with the arbitrary extension to the filename that crab
#now adds, I will add a second associative array to keep track of the
#whole file name. (completefilenames)
#The only subtle part is where the file names are sorted.
#Hopefully this does not create any bugs
indexdict = {}
completefilenames = {}

stub = 'stub'

#do some accounting
nmoved = 0

for line in f:
    mypair = line.split()
#file size
    size = int(mypair[0])
#parse filename
    result = re.match('(.*)_([0-9]*)_([0-9]*)_(.*).root',mypair[1])
    if size > 0:
        if stub=='stub':
            stub = result.group(1)
        elif stub!=result.group(1):
            print "WARNING -- there is more than one group of filenames here!"
        
        if result.group(2) in indexdict:
            #so the key already exists
            indexdict[result.group(2)].append(result.group(3))
            completefilenames[result.group(2)].append(mypair[1])
        else:
            indexdict[result.group(2)] = [result.group(3)]
            completefilenames[result.group(2)] = [mypair[1]]
#        print result.group(0)
#        print result.group(1)
#        print result.group(2)
#        print result.group(3)

#print indexdict

for ii in indexdict:
    if len(indexdict[ii])==1:
        print "nothing to do for ",ii
    else:
        indexdict[ii].sort()
        #will this work? in my one test case, yes
        completefilenames[ii].sort()
        goodindex = indexdict[ii].pop()
        goodfilename = completefilenames[ii].pop()
        print ii, ": keeping index ",goodindex
        print ii, ": corresponds to file ",goodfilename
        for jj in completefilenames[ii]:
            if alreadymadedir == 0:
                print mkdircommand
                #this will give a harmless error if the dir already exists
                os.system(mkdircommand)

            s=inputdir
            s+=jj
            print "moving ", s
            cpcmd = 'rfcp '
            cpcmd += s
            cpcmd += ' '
            cpcmd += outputdir
            print cpcmd
            os.system(cpcmd)
            rmcmd = 'rfrm '
            rmcmd += s
            print rmcmd
            os.system(rmcmd)
            nmoved = nmoved+1

f.close()
os.remove(tmpfile)

print "----------------------------"
print "I moved this many files: ",nmoved
print "----------------------------"
