import os

#Create an executable build script which builds the ma28 library.
#The "gmake clean" ensures we always rebuild the library which is
#necessary to prevent the following example situation: we switch the FLASH
#compiler from absoft to lahey but the library is still built with absoft.
def create_build_script(absLibDir,buildFlag):
    buildScript = absLibDir + '/build.csh'
    if os.path.isfile(buildScript):
      os.remove(buildScript)
    
    fileObj = open(buildScript, 'w')
    fileObj.write('#!/bin/csh\n')
    fileObj.write('#  Dont forget to make this file executable!\n\n')
    fileObj.write('cd source\n')
    fileObj.write('gmake clean\n')
    fileObj.write('gmake BUILDFLAG=' + buildFlag + '\n')
    fileObj.write('cd ..\n')
    fileObj.close()
    os.chmod(buildScript, 0744)


#We use this function to create a custom build file for our internal library.
#This is required because we want to use the compilation flags specified in
#the FFLAGS_[DEBUG|TEST|OPT] variable when building the library. This will
#avoid a prior problem where by the library was not compiled with reals
#promoted to 8 bytes, but the rest of the FLASH source used 8 byte reals.
def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):
    create_build_script(absLibDir,buildFlag)

    #Specify that we want to rebuild the library each time we resetup.
    #Also label the library as internal because we provide the source code.
    #Do not provide a path because we do not have several different
    #ma28 implementations that exist in different directories.
    return {"REBUILD":1, "INTERNAL":""}
