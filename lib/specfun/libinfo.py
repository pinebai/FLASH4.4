import re, os, sys

#Create an executable build script which builds the library
#given by the "lib" argument.
def create_build_script(absLibDir,lib,buildFlag):
    buildScript = absLibDir + '/build.csh'
    if os.path.isfile(buildScript):
      os.remove(buildScript)
    
    fileObj = open(buildScript, 'w')
    fileObj.write('#!/bin/csh\n')
    fileObj.write('#  Dont forget to make this file executable!\n\n')
    fileObj.write('cd ' + lib + '/source\n')
    fileObj.write('make clean\n')
    fileObj.write('make BUILDFLAG=' + buildFlag + '\n')
    fileObj.write('cd ../..\n')
    fileObj.close()
    os.chmod(buildScript, 0744)


def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):
    DEFAULT_LIB = "specfun_subset"
    FULL_LIB = "specfun_dp"

    lib = ""    
    ans = {}
    args = args.lower()
    if not args or re.match("subset", args) or re.match("specfun_subset", args):
       lib = DEFAULT_LIB
    elif re.match("external", args) or re.match("specfun_dp", args):
       lib = FULL_LIB
       return {"EXTERNAL":lib}
    else:
       print >>sys.stderr, 'Unknown specfun library variant "%s"' % args

    if lib:
       create_build_script(absLibDir,lib,buildFlag)
       ans = {"REBUILD":1, "INTERNAL":lib}
    return ans
