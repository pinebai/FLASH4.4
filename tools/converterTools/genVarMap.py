#!/usr/bin/python

import os
import string
import sys
import getopt


usage = '''genVarMap.py [options] <flash2file>
             options:
                -o <filename> : the name of the file that contains the mapping
                                default: flash2vars
                -p [filename] : the name of the file that contains the
                                particle mapping.
                                default name: flash2parts
                -d            : use default values for all mappings
                -h, --help    : prints this message.
                --noPartMap   : do not attempt to generate a particle
                                mapping'''


unkNames = []


try:
    optList, args = getopt.getopt(sys.argv[1:], "o:p:dh",
                                  ["help", "noPartMap"])
except getopt.GetoptError, e:
    print str(e)
    #print "ERROR: Unrecognized option, ", e
    print usage
    sys.exit(1)
except:
    print "unknown error"
    sys.exit(1)
    
#print optList
#print args
doParticles = True
defaultOnly = False
varMapFile = 'flash2vars'
partMapFile = 'flash2parts'

for opt, arg in optList:
    if opt in ('-h', 'help'):
        print usage
        sys.exit(0)
    elif opt == '-o':
        varMapFile = arg
    elif opt == '-p':
        partMapFile = arg
    elif opt == '-d':
        defaultOnly = True
    elif opt == 'noPartMap':
        doParticles = False
    else:
        print "Something has gone very, very wrong. Bailing out."
        sys.exit(1)

#print varMapFile
#print partMapFile
#print defaultOnly
#print doParticles

filename = args[0]

cmd = "h5dump -x " +  filename + " | ./makeVarMap.py"

handle = os.popen(cmd, 'r')
output = handle.readlines()
handle.close()
for line in output :
    unkNames.append(line.strip("\n"))
print unkNames

#outlines now contains a list of our flash2 variable names
activeUnkNames = 0
f3unkNames = []
totVars = 0

if defaultOnly:
    for name in unkNames:
        if name == " 1  ":
            unkNames.remove(name) #this is not included anymore. Originally mfrc
        else:
            f3unkNames.append(name)
            totVars += 1
            
else:  #user wants to check mapping on own
    for name in unkNames :
        if name == " 1  ":
            f3name = raw_input("%s maps to (ommitted by default): " %
                               (name))
        else:
            f3name = raw_input( "%s maps to (default '%s'): " % (name, name))
            #print f3name

        #only add if wanted
        if f3name == "nothing\n": #omitted
            unkNames.remove(name)
        elif f3name == "": #use default
            if name == " 1  ":
                unkNames.remove(name)
                continue
            f3unkNames.append(name)
            totVars += 1
        else:
            f3unkNames.append((f3name + "    ")[0:4])
            totVars += 1

print f3unkNames
f = open(varMapFile, 'w+')
f.write(str(totVars) +"\n")
i = 0
for f2name in unkNames:
    f.write( f2name + " " + f3unkNames[i]+"\n")
    i+=1
f.close()

