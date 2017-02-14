#!/usr/bin/env python
import sys, os, glob, re, string
from Makefile import UnitMakefile
#Add the parent directory to our search path.
sys.path.append(os.path.split(sys.path[0])[0])
from getUseAndCallInfo import getUseAndCallInfo


NON_UNIT_MAKEFILES = [\
"Makefile.h", \
"Makefile.Depend" \
]

#The script parses the Fortran source files in the object directory and
#records every single subroutine / function that is called.  It then
#parses the unit Makefiles in the object directory and records every
#single object file that is needed to build the application binary.
#Now comes the big assumption:
# * Every FLASH subroutine exists in its own file and the file name
#   is the subroutine name plus the .F90 file extension.
#This assumption allows us to remove object files from the Makefile if
#rhe subroutine x in x.F90 (and thus x.o object file) is never
#called in the FLASH application

#There are, of course, sometimes problems making such big assumptions.
#We can add any files that cause problems to IGNORE_FILE_PATTERNS
#list.  Current files to ignore:
# * data - All files containing "data" because most are data modules.
# * interface - All files containing "interface" because most are interfaces.
# * Driver_abortFlashC - never called from a Fortran subroutine, but
# is called from C functions (which we do not parse).
# * RuntimeParameters_[getNum,getPrev,getAll] - overloaded Fortran
# subroutines.
# * gr_pfftFnArg* - Subroutines passed as subroutine arguments.
# * gr_pfftNodeObject, gr_pfftListObject - A linked list data structure.
# * gr_pfftMakePencilIn3dSpace - See gr_pfftGetProcGrid.F90:75
# * gr_hgPoissonSolveBlock - Subroutine passed as subroutine argument.

#Case insensitive.
IGNORE_FILE_PATTERNS = [\
"data", \
"interface", \
"Driver_abortFlashC", \
"RuntimeParameters_getNum", \
"RuntimeParameters_getPrev", \
"RuntimeParameters_getAll", \
"gr_pfftFnArg", \
"gr_pfftNodeObject", \
"gr_pfftListObject", \
"gr_pfftMakePencilIn3dSpace", \
"gr_pfftNodeFnPrototypes", \
"gr_getIndex", \
"gr_packBCs", \
"gr_ptMarkRefineDerefine", \
"gr_hgPoissonSolveBlock" \
]

LOCAL_API_DICT = {\
"gr" : "Grid" \
}


def getUnitNames(objdir):
    nonUnitMakefiles = [objdir+"/"+f for f in NON_UNIT_MAKEFILES]
    allMakefiles = glob.glob(objdir + '/Makefile.*')
    return [x.split(".")[1] for x in allMakefiles if \
            x not in nonUnitMakefiles]


def getNamesOfCalledFns(namePrefixes,sourceFiles):
    calledFns = {}
    for f in sourceFiles:
        print "[getNamesOfCalledFns]: processing file", f
        file = open(f, "r")
        lines = file.read().split("\n")
        file.close()

        for namePrefix in namePrefixes:
            allSIs = getUseAndCallInfo(namePrefix, lines)
            fnList = [x[0] for si in allSIs for x in si.calledSubnames]

            #Store the names of all functions that have a given name prefix.
            #The name prefix is the key and the function name list is the
            #value.  Add previously found functions to the function name list
            #and remove duplicate function names.
            if fnList:
                if calledFns.has_key(namePrefix):
                    fnList.extend(calledFns[namePrefix])
                calledFns[namePrefix] = list(set(fnList))
    return calledFns


#Argumemts:
#logfile - the logfile for writing progress information.
#unitnames, e.g. ["Hydro", "Gravity"] - a list of units that the setup
#                                       script has removed.
#objdir, e.g. PFFT_PoissonFD_compact  - the object directory containing the
#                                       source files which should be modified.
def remove_unit_refs(logfile,unitnames,objdir):

    logfile.write("\n*** Phase 1 ***\n\nEdit Fortran source files by " \
                      "removing references to the following units:\n\n" + \
                      "\n".join(unitnames) + \
                      "\n\nThe removed lines are shown below:\n")
    sWildSearch = objdir + '/*.F90'

    for f in glob.glob(sWildSearch):
        print "[remove_unit_refs]: processing file", f

        file = open(f, "r")
        lines = file.read().split("\n")
        file.close()

        i = 0
        linesToDelete = []

        #Figure out which lines we need to delete from the Fortran
        #source file.  These can be use lines or subroutine calls.
        while i < len(unitnames):
            allSIs = getUseAndCallInfo(unitnames[i], lines)

            for si in allSIs:

                #Handle use lines
                #The indices in si.useOnlyLineIndices is every single line
                #with a used file.  It has taken into account continuation characters
                for a in si.useOnlyLineIndices:
                    linesToDelete.append((a,a))

                #Handle subroutine calls
                for tup in si.calledSubnames:
                    endline = tup[2]
                    while lines[endline].strip().endswith("&"):
                        endline+=1
                    linesToDelete.append((tup[2],endline))

            i+=1

        if (linesToDelete):
            #Sort in reverse order so that we delete from the bottom
            #of the file to the top of the file.  This keeps the line
            #numbers correct.
            linesToDelete.sort(reverse=True)
            logfile.write("\n" + f + "\n")
            for s,e in linesToDelete:
                logfile.write("\n".join(lines[s:e+1]) + "\n")
                del lines[s:e+1]

            lines = "\n".join(lines)
            file = open(f, "w")
            file.writelines(lines)
            file.close()


# We need to fix up the Makefiles so that there are no
# cross-unit dependencies, e.g.
# IO_writeCheckpoint.o : Logfile_stamp.o
# This could happen if a compact application has I/O but no logfiles.
#
def remove_invalid_dependencies(logfile,rmUnits,objdir):

    logfile.write("\n*** Phase 2 ***\n\nEdit Makefiles by " \
                      "removing references to the following units:\n\n" + \
                      "\n".join(rmUnits) + \
                      "\n\nThe invalid dependencies are shown below:\n")

    sWildSearch = objdir + '/*.F90'
    sourceFiles = glob.glob(sWildSearch)
    units = getUnitNames(objdir)
    makeTuples = [(x,objdir + "/Makefile." + x) for x in units]

    for unit,makefileName in makeTuples:
        makefile = UnitMakefile()
        makefile.readFromFile(unit,makefileName)
        makeObj = makefile.getObjList()
        deps = makefile.getDependencies()

        #Obtain a list of foreign-unit object files which appear
        #in the dependency rules.
        depObj = list(set([x for vals in deps.values() for x in vals]))
        unitObjTuple = [(obj.split('_')[0].split('.')[0], obj) for obj in depObj]
        rmDeps = [y for (x,y) in unitObjTuple if x in rmUnits]
        if rmDeps:
            logfile.write("\n" + makefileName + "\n" + \
                              "\n".join(rmDeps) + "\n")
            print "Removing the foreign-unit objects", rmDeps, \
                "from Makefile", makefileName
            makefile.removeObjects(rmDeps)
            makefile.writeToFile(makefileName)


def removeUncalledSubs(logfile,objdir,calledFns):

    logfile.write("\n*** Phase 3 ***\n\nRemove all uncalled Fortran " \
                      "subroutines that have the following name prefix:\n\n" + \
                      "\n".join(unitList) + \
                      "\n\nThe called subroutines are shown below:\n")
    for k in calledFns:
        logfile.write(k + " prefix\n  " + \
                      "\n  ".join(calledFns[k]) + "\n\n")

    #Do not remove a specific object file from the makefile if it
    #matches the pattern in ignorePattern.
    ignorePattern = re.compile("|".join(IGNORE_FILE_PATTERNS), re.IGNORECASE)

    for namePrefix in calledFns:
        if LOCAL_API_DICT.has_key(namePrefix):
            unit = LOCAL_API_DICT[namePrefix]
        else:
            unit = namePrefix
        
        makefileName = objdir + "/Makefile." + unit
        makefile = UnitMakefile()
        makefile.readFromFile(unit,makefileName)
        makeObj = makefile.getObjList()

        #Remove all subroutines from the Makefile which have a given name
        #prefix and are not apparently called (refPattern) and are not
        #ignored (ignorePattern).  The ignore pattern is needed to work
        #around deficiencies of this script and getUseAndCallInfo script.
        unitChars = len(namePrefix)
        refPattern = re.compile("|".join(calledFns[namePrefix]), re.IGNORECASE)
        remList = [x for x in makeObj if \
                   x[0:unitChars] == namePrefix and \
                   not refPattern.search(x) and \
                   not ignorePattern.search(x)]

        #This is just here for debugging:
        logfile.write(namePrefix + " name prefix\n")
        if remList:
            logfile.write(" Will remove the following objects from the " \
                        "Makefile:\n  " + "\n  ".join(remList) + "\n")
            warningObj = [x for x in remList if x not in makeObj]
            if warningObj:
                logfile.write(" Warning! Found objects corresponding to " \
                            "subroutines in the " \
                            "source but not in the Makefile:\n  " +
                            "\n  ".join(warningObj) + "\n")
        logfile.write("\n")


        #Remove the source files and amend the Makefiles.
        sourceList = [objdir + "/" + x.replace(".o",".F90") \
                      for x in remList]
        if sourceList:
            for x in sourceList:
                if os.path.isfile(x):
                    logfile.write(" Removing file " + x + "\n")
                    os.remove(x)
                else:
                    logfile.write(" Cannot remove file " + x + \
                                " because it does not exist\n")
            logfile.write("\n")

        makefile.removeObjects(remList)
        makefile.writeToFile(makefileName)


if __name__=="__main__":
  if len(sys.argv) < 3:
    print "You must supply two arguments:\n\
    1). A file containing all units to remove from the FLASH application\n\
    2). An object directory containing the source files to be removed"
  else:
    fileUnitsToRemove = sys.argv[1]
    objdir = sys.argv[2]

    lines = open(fileUnitsToRemove).read().split("\n")
    units = []
    for line in lines:
        # We treat flashUtilities specially.  All subroutines
        # in flashUtilities have the prefix ut_, but we cannot remove all
        # ut_ because some flashUtilites directories may be needed by the app.
        # We trust the user to remove only those flashUtilities directories
        # that the app does not call.
        if "flashUtilities" not in line:
            h,t = os.path.split(line)
            if t:
                units.append(t)

    print "Removing references to units:", units
    print "The object directory to be edited is:", objdir

    logfile = open(objdir + ".log", 'w')

    #Phase 1
    remove_unit_refs(logfile,units,objdir)

    #Phase 2
    remove_invalid_dependencies(logfile,units,objdir)

    sWildSearch = objdir + '/*.F90'
    sourceFiles = glob.glob(sWildSearch)

    unitList = getUnitNames(objdir)


    #In FLASH a name prefix is either a unit name, e.g. "Grid_"
    #or a localapi name, e.g. "gr_".  We only consider the
    #unit names that are not part of the units selected for removal.
    namePrefixes = [x for x in unitList if x not in units]
    namePrefixes.extend(LOCAL_API_DICT.keys())
        
    #Phase 3
    print "Statically removing unused subroutines that belong to units:", namePrefixes
    publicAPIdict = getNamesOfCalledFns(namePrefixes,sourceFiles)

    #print "The dictionary of called subroutines is:", publicAPIdict
    removeUncalledSubs(logfile,objdir,publicAPIdict)
    logfile.close()
