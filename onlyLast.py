import os, sys, re

#handle input arguments
#this includes the name of the script
##for arg in sys.argv:
##    print arg
    
#sys.argv[1] will have the first argument
tmpfile = '/tmp/joshmt/onlyLast_py_'
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

print mkdircommand
#this will give a harmless error if the dir already exists
os.system(mkdircommand)

f = open(tmpfile,'r')

indexdict = {}

stub = 'stub'

for line in f:
    mypair = line.split()
#file size
    size = int(mypair[0])
#parse filename
    result = re.match('(.*)_([0-9]*)_([0-9]*).root',mypair[1])
    if size > 0:
        if stub=='stub':
            stub = result.group(1)
        elif stub!=result.group(1):
            print "oh no!"
        
        if result.group(2) in indexdict:
            #so the key already exists
            indexdict[result.group(2)].append(result.group(3))
        else:
            indexdict[result.group(2)] = [result.group(3)]
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
        goodindex = indexdict[ii].pop()
        print ii, ": keeping index ",goodindex
        for jj in indexdict[ii]:
            s=inputdir
            s+=stub
            s+='_'
            s+=ii
            s+='_'
            s+=jj
            s+='.root'
#            print "moving ", s
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


f.close()
os.remove(tmpfile)


