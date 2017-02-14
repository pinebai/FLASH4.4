#!/usr/bin/env python

"""This program parses all the .F90 files in the current directory and finds the modules 
each file depends on. Having done that, it generates a makefile listing those dependencies"""

import re, os, os.path, string, sys


usemodpatt = "^\s*use\s+([A-Za-z_0-9]+)[,]?.*"
incpatt = """^\s*[#]\s*include\s+["']?([A-Za-z._0-9]+)['"]?\s+.*"""
defmodpatt = "^\s*module\s+([A-Za-z][A-Za-z_0-9]*)\s*$"
reordpatt = "^\s*!!\s*REORDER[(]([2345])[)]:.*$"
GenmodFilenamepatt = "^.+__genmod$"
#PMdirpatt only works if the source file is a soft-link.
#PMdirpatt = "(PM4_package)|(flash_avoid_orrery)"
PMfilepatt = "^(mpi_|amr_)\w+.F90"

#Code to generate the misc file list (we search all directories included
#by the FLASH Config files):
# for i in $(find {Paramesh4.0,Paramesh4dev}/PM4_package/{headers,source,mpi_source,utilities/multigrid} \
# Paramesh4dev/{flash_avoid_orrery,PM4_package/bittree} -name '*.F90'); do basename $i; done | \
# egrep -v '^(mpi_|amr_).*.F90' | sort | uniq | awk '{ print "\x27" $1 "\x27" }'
PMmiscfiles = ['bittree.F90',
               'clean_divb.F90',
               'clean_field.F90',
               'compress_fetch_list.F90',
               'constants.F90',
               'fill_old_loc.F90',
               'find_surrblks.F90',
               'gr_initParameshArrays.F90',
               'gtest_neigh_data1.F90',
               'io.F90',
               'local_tree_build.F90',
               'local_tree.F90',
               'mesh_test.F90',
               'morton_sort.F90',
               'paramesh_comm_data.F90',
               'paramesh_dimensions.F90',
               'paramesh_interfaces.F90',
               'paramesh_mpi_interfaces.F90',
               'physicaldata.F90',
               'poisson_sor.F90',
               'process_fetch_list.F90',
               'prolong_arrays.F90',
               'quicksort_index.F90',
               'rationalize_fetch_list.F90',
               'rationalize_list.F90',
               'send_block_data.F90',
               'set_f2c_indexes.F90',
               'sparse_solver.F90',
               'test_multigrid.F90',
               'timings.F90',
               'tree.F90',
               'tree_search_for_surrblks.F90',
               'user_coord_transfm.F90',
               'workspace.F90']

# returns the name of a file in current directory with
# almost given name (modulo case changes)
def findMODmatch(name):
    global filelist
    global MODcache
    global MODlist
    if MODlist.has_key(name.lower()):
        MODcache[name] = MODlist[name.lower()] 
    if MODcache.has_key(name): return MODcache[name]
    obj = re.compile(name+"\.[fF](90)?$",re.I)
    for x in filelist:
        m = obj.match(x)
        if m:
           MODcache[name],junk = os.path.splitext(x)
           MODcache[name] = MODcache[name] + ".o"
    if not MODcache.has_key(name): MODcache[name] = ""
    return MODcache[name]

# searches for file "name" in a specified list of directories 
# and returns the absolute name
def findINCmatch(name):
    global INCDirs
    global INCcache
    if INCcache.has_key(name): return INCcache[name]
    for d in INCDirs:
        dname = os.path.join(d,name)
        if os.path.isfile(dname):
           INCcache[name] = dname
           return dname
    INCcache[name]=""
    return ""
 
# Given a filename computes the module and include files it depends on.
# returns a triple of lists:
#       (list of module files  -- a list of strings [not including ".mod"] ,
#        list of include files -- a list of strings [including file extensions] ,
#        list of F90 modules defined in this file -- list of strings [not including ".mod"] ).
def depends(filename):
    """Handle one file"""
    global MODlist, REORDlist, REORDlistPM
    mods = {}
    incs = {}
    heremodules = []
    reorddata = {"FOUR" : [],"FIVE":[],"FLAGS":{} }
    foundREORD = 0
    if not os.path.isfile(filename): return ([],[],[])
    for x in file(filename).readlines():
        m = modRE.match(x)
        if m: 
           mods[m.group(1)] = 1
           continue
        m = incRE.match(x)
        if m: 
           incs[m.group(1)] = 1
           continue
        m = defmodRE.match(x)
        if m: 
           mn = m.group(1)
           basename,ext = os.path.splitext(filename)
           if (mn.lower() == basename) and (basename+".o" == findMODmatch(mn)):
               MODUsed[mn.lower()] = 1
           MODlist[mn.lower()] = basename + ".o"
           heremodules.append(mn.lower())
           continue
        if not foundREORD:
            m = reordRE.match(x)
            if m:
                #paramesh needs special handling.
                if parameshFileRE.search(filename) or filename in PMmiscfiles:
                    REORDlistPM.append(filename)
                else:
                    #not in our list of special directories to account for
                    REORDlist.append(filename)
                foundREORD = 1
                continue
    return ( mods.keys(), [findINCmatch(x) for x in incs.keys()], heremodules )


##################
def main():
   global modRE, incRE, defmodRE, reordRE, parameshFileRE
   global MODcache, INCcache, MODlist, REORDlist, REORDlistPM
   global MODUsed
   global filelist, INCDirs

   # setup global variables
   modRE = re.compile(usemodpatt,re.I)
   incRE = re.compile(incpatt,re.I)
   defmodRE = re.compile(defmodpatt,re.I)
   reordRE = re.compile(reordpatt,re.I)
   parameshFileRE = re.compile(PMfilepatt,re.I)
   genmodFnRE = re.compile(GenmodFilenamepatt,re.I)

   # Process all the files 
   MODcache = {}
   INCcache = {}
   MODlist = {}
   MODUsed = {}
   REORDlist = []
   REORDlistPM = []

   filelist = os.listdir(".") # All files in current directory
   INCDirs = [""] # List of directories which contain include files. "" is current directory
   toProcess = [] # List of files to generate dependency info for
   MODlistObj = {}
   MODlistDone = {}
   MODlistldeps = {}
   MODlistudeps = {}
   MODlistAllHere = {}
   MODlistFortFile = {}
   MODlistCmd = {}

   # Don't generate Makefile lines like '.INTERMEDIATE: Grid_data.mod' unless requested by a flag
   generateINTERMEDIATElines = 0

   # Process argument list and isolate filenames , inc dirs, arbitrary flags 
   for arg in sys.argv[1:]: # for each option given
       if not arg: continue
       if arg[0] == "-" and len(arg) > 2 and arg[1] == "I": # option -I
          INCDirs.append(arg[2:])
       elif arg == "--generateINTERMEDIATElines":
           generateINTERMEDIATElines = 1
       elif arg[0] != "-": # not an option so must be a filename to process
          if arg[0] != '*': toProcess.append(arg) # unprocessed wild card means no such files
       else: 
          pass # Ignore all other options

   ofd = file("Makefile.Depend","w")
   ofd.write("\n###Auto generated file###\n")
   ofd.write("###will be clobbered during the next make\n\n")
   ofd.write("\n# Note that MODUPPERCASE is defined by Makefile.h for sites and compilers where this is necessary,\n")
   ofd.write("# and Makefile should include Makefile.h before Makefile.Depend .\n\n")
   for filename in toProcess:
       basename,ext = os.path.splitext(filename)
       if ext == '.f90' and genmodFnRE.match(basename):
           # Ignoring files *__genmod.f90 which may have been automatically created by Intel compilers
           continue
       mods,incs,modsDefinedHere = depends(filename)
       for h in incs: # Add more files to process if required
           if (h in filelist) and (h not in toProcess): # new header file in current directory
              toProcess.append(h)
       for m in mods:
           MODlist[m.lower()] = findMODmatch(m)
           MODUsed[m.lower()] = 1
       incpart  = " ".join(incs)
       if (not mods) and (not modsDefinedHere): # does not depend on modules 
          ofd.write("%s.o: %s\n" % (basename, incpart))
       else: # depending on MODUPPERCASE use upper or lower case names for .mod files
          if mods and modsDefinedHere: # this file both uses AND defines modules...
              mods = [x.lower() for x in mods] # make sure they are lowercase
              mods = [x for x in mods if (0==modsDefinedHere.count(x))] # strip mods defined here from mods.
          lmodpart = " ".join([x.lower()+".mod" for x in mods])
          umodpart = " ".join([x.upper()+".mod" for x in mods])
          lheremod = " ".join([x.lower()+".mod" for x in modsDefinedHere])
          uheremod = " ".join([x.upper()+".mod" for x in modsDefinedHere])
          if modsDefinedHere:
              firstModDefinedHere = modsDefinedHere[0]
              if (MODlistObj.has_key(firstModDefinedHere)
                  and (MODlistObj[firstModDefinedHere] != firstModDefinedHere+".o")):           # in this case actual writing is deferred to module loop below
                  ofd.write("\nifdef MODUPPERCASE\n%s %s: %s\nelse\n%s %s: %s\nendif\n\n" 
                            % (MODlistObj[firstModDefinedHere], " ".join([x.upper()+".mod" for x in MODlistAllHere[firstModDefinedHere] if x not in modsDefinedHere]), MODlistudeps[firstModDefinedHere],
                               MODlistObj[firstModDefinedHere], " ".join([x.lower()+".mod" for x in MODlistAllHere[firstModDefinedHere] if x not in modsDefinedHere]), MODlistldeps[firstModDefinedHere]))
                  modsDefinedThere = MODlistAllHere[firstModDefinedHere]
                  for m in modsDefinedThere:
                      if MODlistAllHere.has_key(m): del MODlistAllHere[m]
                      if MODlistObj.has_key(m): del MODlistObj[m]
              if (not MODlistObj.has_key(firstModDefinedHere)) or (basename.lower() == firstModDefinedHere):           # in this case actual writing is deferred to module loop below
                  for m in modsDefinedHere:
                      MODlistObj[m] = basename + ".o"
                      MODlistldeps[m] = "%(lmodpart)s %(incpart)s" % locals()
                      MODlistudeps[m] = "%(umodpart)s %(incpart)s" % locals()
                      MODlistFortFile[m] = filename
                      MODlistAllHere[m] = modsDefinedHere
                      if ext == ".F" or ext == ".f":
                          MODlistCmd[m] = "\t$(ECHO-COMPILING)\n\t$(FCOMP) $(FFLAGS) $(FDEFINES) $<\n"
                      elif ext == ".F90":
                          MODlistCmd[m] = "\t$(ECHO-COMPILING)\n\t$(FCOMP) $(FFLAGS) $(F90FLAGS) $(FDEFINES) $<\n"
                      elif ext == ".f90":
                          MODlistCmd[m] = "\t$(ECHO-COMPILING)\n\t$(FCOMP) $(FFLAGS) $(f90FLAGS) $(FDEFINES) $<\n"
              else:
                  ofd.write("\n##ifdef MODUPPERCASE\n##%(basename)s.o %(uheremod)s:%(filename)s %(umodpart)s %(incpart)s \n##else\n##%(basename)s.o %(lheremod)s:%(filename)s %(lmodpart)s %(incpart)s\n##endif\n\n" % locals())
          else:
              ofd.write("\nifdef MODUPPERCASE\n%(basename)s.o %(uheremod)s: %(umodpart)s %(incpart)s \nelse\n%(basename)s.o %(lheremod)s: %(lmodpart)s %(incpart)s\nendif\n\n" % locals())
              

   ofd.write("\n\n### Dependencies of modules\n")
   flash_modules = []
   for (m,obj) in MODlist.items():
       if not MODlistDone.has_key(m):
           if MODlistAllHere.has_key(m):
               if MODUsed.has_key(m): flash_modules.append(MODlistObj[m])
               basename,junk = os.path.splitext(MODlistObj[m])
               allDifferent = 1
               for x in MODlistAllHere[m]:
                   if x.lower() == basename.lower():
                       allDifferent = 0
               ofd.write("\nifdef MODUPPERCASE\n")
               xtraModName = ""
               for x in MODlistAllHere[m]:
                   if (x.upper() != basename): #### and x.upper() == basename.upper()):
                       ofd.write("%s: %s.mod ;if [ -s \"$<\" -a \"$<\" -nt \"$@\" ];then ln -f $< $@;else test -f $@&&touch $@||:;fi\n" %
                                 (x.upper()+".mod",basename))
                       xtraModName = basename + ".mod"
               if xtraModName:
                   if generateINTERMEDIATElines:
                       ofd.write(".INTERMEDIATE: %s.mod\n" % (basename))
                   elif allDifferent:
                       ofd.write(".SECONDARY: %s.mod\n" % (basename))
               ofd.write("%s %s %s: %s %s\n" % (" ".join([x.upper()+".mod" for x in MODlistAllHere[m]]),
                                                xtraModName,
                                                MODlistObj[m],
                                                MODlistFortFile[m],
                                                MODlistudeps[m]) )
               ofd.write("else\n")
               xtraModName = ""
               for x in MODlistAllHere[m]:
                   if (x.lower() != basename): #### and x.lower() == basename.lower()):
                       ofd.write("%s: %s.mod ;if [ -s \"$<\" -a \"$<\" -nt \"$@\" ];then ln -f $< $@;else test -f $@&&touch $@||:;fi\n" %
                                 (x.lower()+".mod",basename))
                       xtraModName = basename + ".mod"
               if xtraModName:
                   if generateINTERMEDIATElines:
                       ofd.write(".INTERMEDIATE: %s.mod\n" % (basename))
                   elif allDifferent:
                       ofd.write(".SECONDARY: %s.mod\n" % (basename))
               ofd.write("%s %s %s: %s %s\n" % (" ".join([x.lower()+".mod" for x in MODlistAllHere[m]]),
                                                xtraModName,
                                                MODlistObj[m],
                                                MODlistFortFile[m],
                                                MODlistldeps[m]) )
               ofd.write("endif\n")
##               if MODlistCmd.has_key(m):
##                   ofd.write(MODlistCmd[m])
               for mdh in MODlistAllHere[m]:
                   MODlistDone[mdh] = 1
           else:
               if MODUsed.has_key(m): flash_modules.append(obj)
               ofd.write("\nifdef MODUPPERCASE\n")
               ofd.write("%s.mod : %s\n" % (m.upper(),obj))
               ofd.write("else\n")
               ofd.write("%s.mod : %s\n" % (m.lower(),obj))
               ofd.write("endif\n")
   ofd.write("\n\n### List of flash module object files\n")
   ofd.write("DATA_OBJS = ")
   for x in range(len(flash_modules)):
       if x % 5 == 0:
          ofd.write("\nDATA_OBJS += ")
       ofd.write(flash_modules[x]+" ")
   ofd.write("\n\n\n")
   ofd.write("%s %s : reorder" % (" ".join(REORDlist), " ".join(REORDlistPM)))
   ofd.write("\n\n\n")

   ofd.close()
   # write out the reorder script
   ifd = open("reorder.tpl")
   ofd = open("reorder.sh","w")
   d = {"reordlist" : " \n ".join(REORDlist), "reordlist_pm" : "\n".join(REORDlistPM)}
   for line in ifd:
       ofd.write(line % d)
   ifd.close()
   ofd.close()

   

if __name__=="__main__":
   main()

