import os, sys, re

###############################################################

#Usage:
# python onlyLast.py <path to input directory>

#now updated to handle output from the latest CRAB with
#file names like this:
#                      myfilename_1_1_aBC.root

# 02/09/11 - Don
# updated to handle files on the T3
# this script now automatically checks if the input dir is pointing
# to the T3 and adapts accordingly
# e.g. python onlyLast.py "$CUSE/blah" 

# 2014 - jmt, remove CASTOR functionality and replace with local
#file commands.

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

#check if the input directory is located on the T3
isT3 = False
checkT3 = 'echo '
checkT3 += inputdir
checkT3 += ' | awk -F "/" \'{print $1}\' '
p = os.popen(checkT3)
location = p.readline()
p.close()
#print location
if "srm" in location:
        print "T3 directory specified.  Using grid commands."
        isT3 = True


outputdir = inputdir
outputdir += 'EXTRAS'

mkdircommand = "mkdir "
mkdircommand += outputdir

#if isT3 == 'T3':
#    mkdircommand = "srmmkdir "
#    mkdircommand += outputdir


#get a file listing into tmpfile
cmd = 'ls -l '
cmd += inputdir
cmd += ' | awk \'/root/ {print $5,$9;}\' > '


if isT3 == True:
    cmd = 'srmls "'
    cmd += inputdir
    cmd += '" | awk -F "/" \'{print $1,$NF}\' | awk \'NF>0\'  > '


cmd+= tmpfile
#print cmd
os.system(cmd)

#keep track of whether we made the directory or not
alreadymadedir = 0
#do not manually make the directory for the T3, lcg-cp will take care of it
if isT3 == True:
    alreadymadedir = 1



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
nmoved_unique = 0
nnotmoved = 0
#check for max file index
max_file_index = -1;

for line in f:
    mypair = line.split()
#file size
    size = int(mypair[0])
#parse filename
    if not mypair[1].endswith('.root'):
        continue #yet another check that the filename ends in .root
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

        #just for user output
        if int(result.group(2)) > max_file_index:
            max_file_index = int(result.group(2))
#        print result.group(0)
#        print result.group(1)
#        print result.group(2)
#        print result.group(3)

#print indexdict

for ii in indexdict:
    if len(indexdict[ii])==1:
        print "nothing to do for ",ii
        nnotmoved +=1
    else:
        nmoved_unique+=1
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
                alreadymadedir=1

            s=inputdir
            s+=jj
            print "moving ", s

            if isT3:
                cpcmd = 'lcg-cp --verbose -b -D srmv2 "'
                cpcmd += s
                cpcmd += '" "';
                cpcmd += outputdir
                cpcmd +='/'
                cpcmd +=jj
                cpcmd += '"'
                print cpcmd
                os.system(cpcmd)
                rmcmd = 'srmrm '
                rmcmd += s
                print rmcmd
                os.system(rmcmd)
            else:
                mvcmd = 'mv %s %s' % (s,outputdir)
                print mvcmd
                os.system(mvcmd)

            nmoved = nmoved+1
    



f.close()
os.remove(tmpfile)

print "----------------------------"
print "Unique files in total: ",nmoved_unique + nnotmoved
print "I moved this many files: ",nmoved
print "I left this many alone: ",nnotmoved
print "Max file index: ",max_file_index
print "----------------------------"



