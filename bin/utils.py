
# Bunch of small functions doing useful stuff

__all__ = [ "updateAndMergeVariablePropertyTuples", "getRelPath", "dirGlob", "stripComments", "determineMachine", 
            "is_upper",  "strictlyCaseSensitiveFilenames",
          ] 

##############################################

import globals
from globals import *  # GVars and SetupError
from lazyFile import * # for LazyFile

import sys,os.path, string, glob, socket, re

def is_upper(letter):
    return letter == letter.upper()

def updateAndMergeVariablePropertyTuples(dictout,dictin,infotext):
    for (k,v) in dictin.items():
        if not dictout.has_key(k):
            dictout[k] = v
        elif type(v)==type([]) or type(v)==type(()):
            def checkIfBothSignificantAndDifferent(l):
                if not l[1] or l[1]=='NONEXISTENT' or l[1]=='GENERIC':
                    return False
                if not l[0] or l[0]=='NONEXISTENT' or l[0]=='GENERIC':
                    return False
                return (l[0] != l[1])

            if not (filter(lambda x:x=='NONEXISTENT' or x=='GENERIC', v)):
##                print "Merging", dictout[k], "<=", v 
                zippedItem = zip(dictout[k], v)
                try:
                    firstProblem = map(checkIfBothSignificantAndDifferent, zippedItem).index(True)
                except ValueError:
                    firstProblem = -1
                if (firstProblem >= 0): # Put some effort into generating meaningful messages.
                    if infotext=="VARIABLE":
                        if firstProblem==0:
                            attr = "TYPE"
                        elif firstProblem==1:
                            if zippedItem[1][0]==zippedItem[2][0] and zippedItem[1][1]==zippedItem[2][1]:
                                attr = "EOSMAP"
                            else:
                                attr = "EOSMAPIN"
                        elif firstProblem==2:
                            attr = "EOSMAPOUT"
                        else:
                            attr = "EOSMAP"
                    else:
                        if firstProblem==0:
                            if zippedItem[0][0]==zippedItem[1][0] and zippedItem[0][1]==zippedItem[1][1]:
                                attr = "EOSMAP"
                            else:
                                attr = "EOSMAPIN"
                        elif firstProblem==1:
                            attr = "EOSMAPOUT"
                        else:
                            attr = "EOSMAP"
                    raise SetupError('Conflicting specifications for %s %s %s: "%s" and "%s".' %
                                     (infotext,k,attr, zippedItem[firstProblem][0],zippedItem[firstProblem][1]))
                dictout[k] = v
            else:
                def secondIfSignificant(l):
                    if not l[1] or l[1]=='NONEXISTENT' or l[1]=='GENERIC':
                        if l[0] and l[0]!='NONEXISTENT' and l[0]!='GENERIC':
                            return l[0]
                    return l[1]
                    
##                print "Merging", dictout[k], "<-", v 
                zippedItem = zip(dictout[k], v)
                try:
                    firstProblem = map(checkIfBothSignificantAndDifferent, zippedItem).index(True)
                except ValueError:
                    firstProblem = -1
                if (firstProblem >= 0): # Put some effort into generating meaningful messages.
                    if infotext=="VARIABLE":
                        if firstProblem==0:
                            attr = "TYPE"
                        elif firstProblem==1:
                            if zippedItem[1][0]==zippedItem[2][0] and zippedItem[1][1]==zippedItem[2][1]:
                                attr = "EOSMAP"
                            else:
                                attr = "EOSMAPIN"
                        elif firstProblem==2:
                            attr = "EOSMAPOUT"
                        else:
                            attr = "EOSMAP"
                    else:
                        if firstProblem==0:
                            if zippedItem[0][0]==zippedItem[1][0] and zippedItem[0][1]==zippedItem[1][1]:
                                attr = "EOSMAP"
                            else:
                                attr = "EOSMAPIN"
                        elif firstProblem==1:
                            attr = "EOSMAPOUT"
                        else:
                            attr = "EOSMAP"
                    raise SetupError('Conflicting specifications for %s %s %s: "%s" and "%s".' %
                                     (infotext,k,attr, zippedItem[firstProblem][0],zippedItem[firstProblem][1]))
                mergedItem = map(secondIfSignificant, zippedItem)
##                print "Mergeditem is", mergedItem 
                dictout[k] = mergedItem
        else:
            dictout[k] = v
                    

def getRelPath(filename,basedir):
    """Return the relative path to FILENAME from basedir"""
    sep = os.sep
    srcDir = os.path.abspath(basedir)+sep
    tgtDir = os.path.abspath(os.path.dirname(filename))+sep
    # most common prefix (will end with "/" or contain extra chars
    cp = os.path.dirname(os.path.commonprefix([srcDir,tgtDir]))
    # Handles the special case when common prefix is "/"
    if cp[-1] != sep: cp = cp + sep
    src = srcDir[len(cp):]
    tgt = tgtDir[len(cp):]
    # src and tgt contains name relative to cp
    c = src.count(sep) # how many levels up to reach common dir
    if c == 0:
      prefix = "."+sep
    else: prefix = c*(".."+sep)
    return os.path.join(prefix+tgt,os.path.basename(filename))


# Takes a pattern (absolute or relative to current directory)
# and returns a list of directories matching pattern. The match is 
# made case-insensitive
def dirGlob(pathname):
    # maps a -> [aA] but "/" -> "/", "1" -> "1" 
    mapfn = lambda x: (x.upper() != x.lower() and "[%s%s]" % (x.lower(),x.upper())) or x
    globstr = string.join(map(mapfn,pathname),"") # concatenate
    # find files which match pathname (except for case)
    files = glob.glob(globstr)
    # return only those of which are directories
    return [name for name in files if os.path.isdir(name)]  

# return the string after removing comments
def stripComments(line, commentChar='#',quoteChar='"'):
    cut = 0
    while cut < len(line):
       pos = string.find(line[cut:],commentChar)
       # no real comment in line
       if pos < 0: 
          return line  
       # do we have even number quotes before it? If so this is real comment
       if string.count(line[:cut+pos],quoteChar) % 2 == 0:
          return line[:cut+pos]
       cut = cut + pos+1
    return line

def getOSType(prototypesDir):
    ostype = string.lower(sys.platform)
    if '-' in ostype:
        ostype = ostype[:string.find(ostype, '-')]
    for proto in os.listdir(prototypesDir):
        if string.count(ostype, string.lower(proto)):
            return proto
    return ostype

# return the hostname to use
def getHostName(sitesDir):

    #Change made by sam 09/01/2011. A bug occurs when trying to get the 
    #host by address for machines without dns entries in a name server 
    #or on a machine without a connection to a dns server
    
    tempHostName = socket.gethostname()

    try:	
    	temp = socket.gethostbyaddr(tempHostName)
    except:
	temp = (socket.getfqdn(),[])
        
	
    fallback = temp[0]
    namesToTry = [temp[0]]
    namesToTry.append(socket.gethostname())
    namesToTry.extend(temp[1]) # list of addl names for the current host

    # Read the alias file into memory
    aliasLines = []
    try:
       aliasFile = open(os.path.join(sitesDir,'Aliases'))
       GVars.out.put('checking sites Aliases file',globals.IMPINFO)
       for line in aliasFile.readlines():
           line = stripComments(line, '#','"')
           line = string.strip(line)
           if line:
              parts = string.split(line)
              if len(parts) <> 2:
                 GVars.out.put("Ignoring bad Aliases file line '%s'" % string.strip(line),globals.WARN)
              else: aliasLines.append(parts)
       aliasFile.close()
    except IOError:
        GVars.out.push()
        GVars.out.put("couldn't open sites Aliases file",globals.WARN)
        GVars.out.pop()

    # try all these hostnames and return the first one which 
    # succeeds. If all fails just return what we would normally
    # have returned
    for name in namesToTry:
        touse = getHostNameToUse(sitesDir,name,aliasLines)
        if touse: return touse
    return fallback

def getHostNameToUse(sitesDir,hostname,aliasLines):
    
    for (site,regex) in aliasLines:
        if re.match(regex, hostname) != None:
           hostname = site
           break
            
    ans = None
    for site in os.listdir(sitesDir):
        if string.count(hostname, site):
            ans = site
    return ans

def determineMachine():
    """Returns directory of proper machine to use"""
    GVars.out.put('checking for needed files and directories',globals.IMPINFO)
    GVars.out.push()
    
    siteDir = os.path.join(GVars.flashHomeDir, 'sites')
    systemsDir = os.path.join(siteDir, 'Prototypes')
    ostype = getOSType(systemsDir)

    if GVars.build_tau:
        GVars.out.put('using TAU stub makefile '+GVars.build_tau,globals.IMPINFO)
    
    if GVars.build_site:
        machDir = os.path.join(siteDir, GVars.build_site)
        if os.path.isdir(machDir):
            GVars.out.put('using site directory for site '+GVars.build_site,globals.IMPINFO)
                          
        else:
            raise SetupError('fatal:  could not find site directory for '\
                             'site %s'%GVars.build_site,globals.ERROR)
    
    elif GVars.build_os:
        machDir = os.path.join(systemsDir, GVars.build_os)
        if os.path.isdir(machDir):
            GVars.out.put('using prototype directory for ostype '+ \
                          GVars.build_os, globals.IMPINFO)
                          
        else:
            raise SetupError('fatal:  could not find prototype directory for '
                             'ostype ' + GVars.build_os, globals.ERROR)

    else:
        hostname = getHostName(siteDir)
        machDir = os.path.join(siteDir, hostname)
        if os.path.isdir(machDir):
            GVars.out.put('using site directory for site '+hostname, globals.IMPINFO)
        else:
            machDir = os.path.join(systemsDir, ostype)
            if os.path.isdir(machDir):
                GVars.out.put('site directory for site '+hostname+\
                              ' not found;',globals.WARN)
                GVars.out.put('using prototype '+ostype,globals.WARN)
            else:
                raise SetupError('fatal:  could not find site for prototype'\
                                 ' directory!\n'
                                 '         specify site or ostype, or else '\
                                 'create a directory for your site\n\n'
                                 '         site    = %s'
                                 '\n         ostype  = %s'%(hostname,ostype))

    GVars.out.pop()
    return machDir

def strictlyCaseSensitiveFilenames():
    """Determines whether case is strictly significant in filenames"""
    GVars.out.put('checking case-sensitivity of filenames',globals.INFO)
    GVars.out.push()

    testname1 = 'TestFileName_tempFile123.mod'
    testname2 = 'testfilename_tempfile123.mod'

    # Strategy:
    # (1) Make sure file testname2 does nto exist in the object directory
    # (2) Generate file testname1 in the object directory
    # (3) If now file testname2 exists in the object directory, filenames
    #     are not handled in a strictly case-independenty manner; otherwise,
    #     assume that they are.
    file1 = os.path.abspath(testname1)
    file2 = os.path.abspath(testname2)

    if os.path.exists(file2):
        os.remove(file2)
    open(file1,'w')
    if os.path.exists(file2):
        ans = 0
    else:
        ans = 1

    os.remove(file1)

    GVars.out.put('determined case sensitivity to be %d' % ans, globals.DEBUG)
    GVars.out.pop()
    return ans

