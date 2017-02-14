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
    fileObj.write('gmake clean\n')
    fileObj.write('gmake BUILDFLAG=' + buildFlag + '\n')
    fileObj.write('cd ../..\n')
    fileObj.close()
    os.chmod(buildScript, 0744)


def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):
    DEFAULT_LIB = "default_sqrt3"
    BG_APPROX_LIB = "bg_approx_sqrt3"

    lib = ""    
    ans = {}
    args = args.lower()
    if not args or re.match("default", args):
       lib = DEFAULT_LIB
    elif re.match("bg_approx", args):
       lib = BG_APPROX_LIB
    else:
       print >>sys.stderr, 'Unknown cube root variant "%s"' % args

    if lib:
       create_build_script(absLibDir,lib,buildFlag)
       ans = {"REBUILD":1, "INTERNAL":lib}
    return ans
