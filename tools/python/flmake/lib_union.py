"""Code for LibUnion class"""

import string
import re
import os 
import UserDict
import sys
import shutil

from . import setup_globals
from .setup_globals import gvars, SetupError
from .utils import get_rel_path
from .flash_lib import FlashLib
from .lazy_file import LazyFile 


class LibUnion(UserDict.UserDict):

   def __init__(self,libs): # list of lib names
       UserDict.UserDict.__init__(self)
       self.adjust(libs)

   def adjust(self,liblist): # satisfy dependencies and perform checks
       pwd = os.getcwd()
       os.chdir(gvars.libDir)
       ans = self.libDepHelper(liblist,[]) 
       # now self["libname"] -> FlashLib("libname")
       os.chdir(pwd)
       self.topoSort(ans)
       # now self.libOrder is a list of pairs containing libname and arguments

   def libDepHelper(self,libs,done):
      """Add libraries which those in libs depend on.
      done is a list of names of libraries already processed"""
      tocheck = None
      for a in libs.keys():
        if a not in done: 
           tocheck = a
           break
      # Nothing left to do
      if not tocheck: return libs
      lib = tocheck
      self[lib] = FlashLib(lib) # get lib Config Info. Handles, internal and external libraries correctly
      newdeps = self[lib]['LIBRARY']
      for x in newdeps.keys(): # for each dependency update libs with args
        if not libs.has_key(x): libs[x] = []
        libs[x].append(newdeps[x])
      done.append(tocheck)
      return self.libDepHelper(libs,done) # recursive call to handle the rest

   # libs is a dictionary mapping libs -> list of arglists (one for each person depending on this lib)
   # self is also a dictionary mapping library name to FlashLib instance for specific library
   #   in particular self["libname"]["LIBRARY"] contains other libs it depends on
   # returns all libaries in an order so that all the dependencies of lib comes before the lib itself
   def topoSort(self,libs):
      curr = [] # list of libraries whose dependencies have been satisfied
      rest = libs.keys() # those whose dependencies need to be satisfied
      while rest: # As long as their is something to do
            success = 0
            for lib in rest:
                libOK = 1 # have we satisfied all dependencies of "lib"
                # self[lib][LIBRARY] maps libname to argument for lib
                for other in self[lib]['LIBRARY'].keys():
                    if other not in curr: libOK = 0
                if libOK == 1:
                   success = 1
                   curr.append(lib)
                   rest.remove(lib)
            if success == 0: # could satisfy any more dependencies
               break
      if rest: # dependency cycle detected
        gvars.out.put("Topological Sort Failed: Cycle Detected",setup_globals.ERROR)
        gvars.out.put("Suspect Libraries: %s"% string.join(rest,","),setup_globals.ERROR)
        raise SetupError("Check your library dependencies")
      # Now check if different people depend on the library in a different way
      # if so flag an error
      self.libOrder = []
      for lib in curr:
          args = filter(None,libs[lib]) # pull out all non-trivial arguments
          if not args:
             arg = ""
          else:
             arg = args[0]
             if [x for x in args if x != arg]: # if there is any entry in args other than args[0]
                gvars.out.put("non-unique non-empty argument lists for library %s\n"%lib,setup_globals.ERROR)
                gvars.out.put("Argument Lists are: %s\n"% str(args),setup_globals.DEBUG)
                raise SetupError()
          self.libOrder.append((lib,arg))

   ######################### remaining methods deal with generating flags for libraries 
   ######################## and building the library if reqd

   # set the list of macros found in Makefile.h
   def setMacros(self,macros=[]):
       self.macros = macros


   def getLibFlags(self,lib,buildFlag,args="",makefilename=""):
     """Returns dictionary of flags to be added to the appropriate stage of building 
     the executable. 

     macros is the list of macros defined in Makefile.h
     args is the list argument strings (ideally there will 0/1 only)

       CFLAGS -> Include flags
       FFLAGS -> include flags
       LIB    -> Linking flags
 
     or if lib is not internal calls extLibFlags to find the flags

     Does not handle the libraries that this library depends on.
     Raises errors if internal lib exists but is troublesome"""


     base = string.lower(lib)
     libDir = os.path.join(gvars.libDir,base)
     #relLibDir = get_rel_path(libDir,".")
     relLibDir = get_rel_path(libDir, gvars.project_setup_dir)
     if not os.path.isdir(libDir): 
        return self.extLibFlags(lib,buildFlag) # args are of no use here

     # does ...lib/libname/libinfo.py exist?
     if os.path.isfile(os.path.join(libDir,'libinfo.py')):
       # use the new method of getting info
       sys.path.insert(1,libDir) # ask python to search here
       # call the libinfo function in libinfo.py and store its result
       libFlags = __import__("libinfo").libinfo(absLibDir=libDir,
                              relLibDir=relLibDir,
                              buildFlag=buildFlag,
                              args=args,
                              macros=self.macros)
       sys.path.remove(libDir)
       del sys.modules["libinfo"]
       # if we are an internal library and need to rebuild
       if libFlags == None: # programming error
          raise SetupError("libinfo for %s returned nothing. Programming error" % base)
       if libFlags.has_key("INTERNAL"): # pick up default flag info
          libFlags = self.intLibFlags(base,subdir=libFlags["INTERNAL"],absLibDir=libDir,relLibDir=relLibDir)
       if libFlags.has_key("EXTERNAL"): # pick up default flag info for external
          libFlags = self.extLibFlags(libFlags["EXTERNAL"],buildFlag)
       if self[lib]["TYPE"]=="INTERNAL" and libFlags.get("REBUILD",None):
          self.makeBinary(libDir,base,makefilename)
       return libFlags

     # now we know that there is no libinfo.py
     if self[lib]["TYPE"]=="EXTERNAL": 
        return self.extLibFlags(lib.upper(),buildFlag)

     # no libinfo and internal: Compile the binary if reqd
     if self.libraryExists == False:
        self.makeBinary(libDir,base,makefilename)
        
     # return the default infor for simple internal binaries
     return self.intLibFlags(base,absLibDir=libDir,relLibDir=relLibDir)

   # Return the usual flags
   # libname = name of the library e.g. pfft
   # subdir = subdirectory containing what we want
   # libDir = base directory for library "FLASHHOME/lib/pfft" 
   def intLibFlags(self,libname,subdir="",absLibDir="",relLibDir=""):
    ans = {}
    if subdir != "" :
       relLibDir = os.path.join(relLibDir,subdir)
       absLibDir = os.path.join(absLibDir,subdir)
       libname = subdir.lower()
    else:
       libname = libname.lower()

    ans["LIB"] = '-L%s/object/ -l%s'%(relLibDir,libname)
    if os.path.isdir(os.path.join(absLibDir, 'include')):
        includeMacro = '-I%s/include'%relLibDir
    else:
        includeMacro = ''
    ans["CFLAGS"] = includeMacro
    ans["FFLAGS"] = includeMacro
    if not os.path.isfile(os.path.join(absLibDir,"lib%s.a"%libname)):
       ans["REBUILD"]=1
    return ans

   # return CFLAGS, FFLAGS,... from the makefile
   def extLibFlags(self,lib,buildFlag):
     ans = {}
     for compiler in setup_globals.COMPILERS:
        for macro in ['%s_%s_%s'%(compiler, lib.upper(), buildFlag),
                      '%s_%s_%s'%(compiler, lib.upper(), setup_globals.DEFLTFLAG),
                      '%s_%s'%(compiler, lib.upper())]:
            if macro in self.macros:
               ans[compiler] = '$(%s)'%macro
               break
     return ans

   def libraryExists(self,libDir,libname):
     obj_dir = os.path.join(libDir, 'object')
     binary = os.path.join(obj_dir, 'lib%s.a'%libname)
     if os.path.isfile(binary):
        return True
     else:
        return False

   # Build the binary if required
   def makeBinary(self,libDir,libname,makefilename):
     USAGE = """
     Fatal Error: setup is unable to locate the internal library %(libname)s.
     In order to have setup identify an internal library, one of the following should hold:

     * the file %(libDir)s/libinfo.py must exist and contain a function libinfo
       this function must return a dictionary containing flags to be passed to the
       different compilers and linkers. libinfo can request that the binary be
       rebuilt by setting the REBUILD flag in its return value
     * the file %(libDir)s/object/lib%(libname)s.a must exist 
     * the file %(libDir)s/build.csh or %(libDir)s/build.py must exist
       this script is given the location of Makefile.h as its only argument
       - If there is a libinfo.py the build script must generate 
         the library specified by return value of libinfo function
       - Otherwise the build script must generate %(libDir)s/object/lib%(libname)s.a
     """
     pwd = os.getcwd()
     os.chdir(libDir)
     # Link/Copy current Makefile.h so library is built using
     # the same compiler and options as the FLASH code
     if os.path.isfile("Makefile.h"): os.remove("Makefile.h")
     if gvars.portable:
        shutil.copy2(makefilename,"Makefile.h")
     else:
        os.symlink(makefilename,"Makefile.h")
        
     if os.path.isfile('build.py'):
        gvars.out.put('exec\'ing %s/build.py'%libDir)
        exec(open('./build.py').read())
     elif os.path.isfile('build.csh'):
        gvars.out.put('running %s/build.csh'%libDir)
        os.system('./build.csh')
     else:
        raise SetupError(USAGE % locals())
     os.chdir(pwd)

   def write_libraries(self):
    #fd = LazyFile(os.path.join(gvars.flash_src_dir,gvars.project_build_dir,setup_globals.SETUP_LIBRARIES_FILENAME))
    fd = LazyFile(os.path.join(gvars.project_setup_dir, setup_globals.SETUP_LIBRARIES_FILENAME))
    for (libname,args) in self.libOrder:
        if args:
           fd.write("%8s : %s\n" % (libname,args))
        else: fd.write("%8s : \n" % libname)
    fd.close()
    if not fd.samefile and gvars.noClobber:
       gvars.out.put("WARNING: Library requirements have changed from previous run. Ignoring noClobber")
       gvars.noClobber = 0


