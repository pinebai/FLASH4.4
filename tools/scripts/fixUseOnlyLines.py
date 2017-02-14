#!/usr/bin/env python
import sys, os
from getUseAndCallInfo import getUseAndCallInfo
from interfacesDict import *
delStr = ">>>>>DELETE THIS LINE<<<<<"

def adjust(files, interfacesDict, unitnames, dryRun):

  changedFilesDict = {}
  ctr              = 0
  numFiles         = len(files)

  # cull out any unit-names that don't have an API-level interface
  for unitname in unitnames:
    if interfacesDict.getPathAndSubroutinesForUnit(unitname)[0] == None:
      unitnames.remove(unitname)
      if dryRun:
        print "Error: No API-level interface found for unit \"%s\" - removing this unit from consideration." % unitname

  if not unitnames:
    print "Error: No units to check."
    sys.exit(1)

  # else
  for f in files:
    lines              = open(f).read().split("\n")
    rewrite            = False
    linesNeedAdjusting = False
    ctr += 1

    #CD (18 Apr 08): This progress indicator is always useful, whether it be dryRun or not.
    print "working on file %s/%s" % (ctr,numFiles)

    # Use "while" instead of typical "for" since there can be
    # use-cases where we examine the same file multiple times.
    # In that case the iterator 'i' is not advanced.
    i = 0
    while i < len(unitnames):
      allSIs = getUseAndCallInfo(unitnames[i], lines)

      for si in allSIs:
        subname               = si.subname
        subroutineLineIndex   = si.subroutineLineIndex
        poundIncludeLineIndex = si.poundIncludeLineIndex
        useOnlyLineIndices    = si.useOnlyLineIndices
        otherUseLineIndex     = si.otherUseLineIndex
        indent                = si.indent
        useOnlySubnames       = si.useOnlySubnames
        calledSubnames        = si.calledSubnames 

        if poundIncludeLineIndex >= 0:
          linesNeedAdjusting = True

        if len(useOnlyLineIndices) > 0 and len(useOnlySubnames) == 0:
          # a 'use interface' is here even though
          # there are no subroutines following it
          linesNeedAdjusting = True

        ## CAPITALIZATION

        # If this subroutine is itself an API-level subroutine, make
        # sure its capitalization is correct
        correctlyCasedMainSubname = interfacesDict.getCorrectlyCasedSubname(subname)
        if correctlyCasedMainSubname and correctlyCasedMainSubname != subname:
          linesNeedAdjusting

        # Examine all subroutines references in the "USE, ONLY" line.
        # If they are API-level subroutines, make sure capitalization
        # is correct. If they aren't, they shouldn't be here.
        correctlyCasedUseOnlySubnames = []
        for useOnlySubname in useOnlySubnames:
          correctlyCasedUseOnlySubname = interfacesDict.getCorrectlyCasedSubname(useOnlySubname)
          if correctlyCasedUseOnlySubname:
            correctlyCasedUseOnlySubnames.append(correctlyCasedUseOnlySubname)
            if correctlyCasedUseOnlySubname != useOnlySubname:
              linesNeedAdjusting = True
          else:
            linesNeedAdjusting = True

        # Examine all called subroutines in the body of the text. If
        # they are API-level subroutines, make sure capitalization is
        # correct. If they aren't, disregard them. They are probably
        # internal subroutines that are case-indistinguishable from
        # API-level ones.
        badCalledSubnames = []
        correctlyCasedCalledSubnames = []
        for tup in calledSubnames:
          calledSubname = tup[0]
          correctlyCasedCalledSubname = interfacesDict.getCorrectlyCasedSubname(calledSubname)
          if correctlyCasedCalledSubname:
            correctlyCasedCalledSubnames.append(correctlyCasedCalledSubname)
            if correctlyCasedCalledSubname != calledSubname:
              badCalledSubnames.append(tup)
              linesNeedAdjusting = True


        # CONSISTENCY

        #------------------------------------------------------------------
        #START
        #CD (18 Apr 08): The removal of names from the use list does not
        #work reliably.  There is no problem with the commented out code below, rather it
        #is a problem with the information gathering phase.  It should be noted
        #that it gathers MOST of the available information, which allows it
        #to add the majority of missing API subroutine names to the use list.
        #However, without complete information it is too risky to allow the
        #script to remove names.
        #   *** It is not perfect, but it is a VERY helpful script ***
        #------------------------------------------------------------------
        
        # Make sure eveything in the "USE, ONLY" line
        # is actually called in the body of the file.
        #for correctlyCasedUseOnlySubname in correctlyCasedUseOnlySubnames:
        #  if correctlyCasedUseOnlySubname not in correctlyCasedCalledSubnames:
        #    correctlyCasedUseOnlySubnames.remove(correctlyCasedUseOnlySubname)
        #    linesNeedAdjusting = True
        #END
        #------------------------------------------------------------------


        # Make sure all called subroutines in the body
        # are present in the "USE, ONLY" line.
        for correctlyCasedCalledSubname in correctlyCasedCalledSubnames:
          if correctlyCasedCalledSubname == correctlyCasedMainSubname:
            continue  # don't include this recursive subroutine's name in its own "USE, ONLY" line
          if correctlyCasedCalledSubname not in correctlyCasedUseOnlySubnames:
            correctlyCasedUseOnlySubnames.append(correctlyCasedCalledSubname)
            linesNeedAdjusting = True


        if linesNeedAdjusting:
          rewrite = True

          # Note which file, subroutine, and unit we're changing
          # for later reporting.
          if changedFilesDict.has_key(f):
            if changedFilesDict[f].has_key(subname):
              changedFilesDict[f][subname].append(unitnames[i])
            else:
              changedFilesDict[f][subname] = [unitnames[i]]
          else:
            changedFilesDict[f] = {subname: [unitnames[i]]}

          if correctlyCasedMainSubname and correctlyCasedMainSubname != subname:
            lines[subroutineLineIndex] = lines[subroutineLineIndex].replace(subname, correctlyCasedMainSubname)

          if poundIncludeLineIndex >= 0:
            lines[poundIncludeLineIndex] = delStr

          for j in useOnlyLineIndices:
            lines[j] = delStr

          # Correct any mis-capitalizations in API-level function calls
          # in the body of the file.
          for badCalledSubname, badCalledSubnameInContext, badLine in badCalledSubnames:
            correctedSubname = interfacesDict.getCorrectlyCasedSubname(badCalledSubname)
            correctedContext = badCalledSubnameInContext.replace(badCalledSubname, correctedSubname)
            lines[badLine] = lines[badLine].replace(badCalledSubnameInContext, correctedContext)

          if len(correctlyCasedUseOnlySubnames) > 0:
            newLine = (" "*indent) + ("use %s_interface, ONLY : " % unitnames[i]) + ", ".join(correctlyCasedUseOnlySubnames)
            if len(newLine) > 70:
              startMark = 0
              while (len(newLine) - startMark > 70):
                lineBreak = newLine.rfind(" ",startMark,startMark+70)+1
                newLine = newLine[:lineBreak] + "&\n" + (" "*(indent+2)) + newLine[lineBreak:]
                startMark += (lineBreak+2)

            if len(useOnlyLineIndices) > 0:
              # this file had a 'use Grid_interface' before, so we put
              # the new line back in where the old one started.
              # Note that 'newLine' might contain '\n' characters, but
              # this does not yet add to the length of 'lines', as we
              # insert all of 'newLine' into a single slot in 'lines'
              lines[useOnlyLineIndices[0]] = newLine

            # up to this point we have carefully *NOT* changed the length
            # of 'lines' by replacing lines to be deleted with 'delStr'
            # and inserting 'newLine' into a single slot in 'lines'.
            # Had we done this, it would have thrown off the various index
            # values we got back from "getUseAndCallInfo()"
            # Now, however, we do actually insert a line. This should be
            # OK, since we've already marked all our lines to be deleted
            # with 'delStr'
            elif otherUseLineIndex >= 0:
              # this file didn't have a 'use Grid_interface' before,
              # but it did have some other 'use' statement, so we slip
              # our new line right in under that
              lines.insert(otherUseLineIndex, newLine)

            elif subroutineLineIndex >= 0:
              # this file didn't have any use statements, so add the
              # new line right under the 'subroutine' declaration
              lines.insert(subroutineLineIndex, newLine)

          # shake out multiple '\n's that might exist inside a single line
          # (such as those that might be introduced when writing 'newLine')
          lines = "\n".join(lines).split("\n")

          # go through once more to delete instances of 'delStr'
          j = 0
          while j < len(lines):
            if lines[j] == delStr:
              del lines[j]
            else:
              j += 1

          # It is possible that we have now changed the number of lines in the file.
          # If the file has more than one subroutine in it, the line numbers that we
          # got back from "getUseAndCallInfo.py" will no longer be accurate below the
          # point at which we made our changes. Therefore we re-submit 'lines' to the
          # script to ensure accurate numbers for other subroutines, if there are any.
          linesNeedAdjusting = False
          break  # 'i' will not advance
      else:
        # all subroutines are now correct
        i+=1
    # out of the inner loop (the while-loop)
    if rewrite:
      if dryRun:
        print "---------" + ("-"*len(f))
        print "file is: %s" % f
        print "---------" + ("-"*len(f))
        print "\n".join(lines)
        print ""
      else:
        open(f,"w").write("\n".join(lines))

  # out of the outer loop (the for-loop)
  if len(changedFilesDict) > 0:
    reportLines = []
    for changedFile in changedFilesDict.keys():
      reportLines.append("file: %s" % changedFile)
      subnamesDict = changedFilesDict[changedFile]
      for subname in subnamesDict.keys():
        reportLines.append("  subroutine: %s" % subname)
        unitnames = subnamesDict[subname]
        reportLines.append("    affected unit(s): " + ", ".join(unitnames))

    # summary of all activity
    if dryRun:
      print "The following files require changes:"
      print "------------------------------------"
    else:
      print "The following files were changed:"
      print "---------------------------------"
    for reportLine in reportLines:
      print reportLine
  else:
    if dryRun:
      print "No files require changes."
    else:
      print "No files were changed."

def usage():
  print "./fixUseOnlyLines.py [-dry-run] [-auto | <unit-name1> <unit-name2>]..."
  print ""
  print "<unit-name> should be the name of a FLASH unit (e.g. Grid). If"
  print "none are provided, the script will check all units for which a"
  print "corresponding _interface.F90 file exists."
  print ""
  print "This script ensures that all files that make calls to subroutines"
  print "beginning with \"<unit-name>_\" declare those subroutines in a"
  print "\"use\" statement below their subroutine header of the form:"
  print ""
  print "  use <unit-name>_interface.F90, ONLY: <call1>, <call2>, etc."
  print ""
  print "It also ensures that any subroutines references in an already-"
  print "existing \"use\" statement are indeed called in the body of the"
  print "code, and removes them from the \"use\" statement if they are not."
  print ""
  print "The script also seeks out and removes instances of:"
  print ""
  print "  #include <unit-name>_interface.h"
  print ""
  print "which is a deprecated style."
  print ""
  print "When run with the \"-dry-run\" flag, the script prints the lines"
  print "of any adjusted files to stdout, but does not actually overwrite"
  print "the file."
  sys.exit(0)

if __name__ == "__main__":
  args = sys.argv
  dryRun = False

  if len(args) == 1:
    usage()

  for item in ["-h", "--h", "-help", "--help"]:
    if item in args:
      usage()

  # else
  for item in ["--dry-run", "-dry-run", "--dry", "-dry", "-d", "--d"]:
    if item in args:
      dryRun = True
      args.remove(item)
      break

  # This script assumes it is being run from "FLASH3/tools/scripts"
  pathToFlash = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(args[0]))))
  f90Files = []
  for root, dirs, files in os.walk(os.path.join(pathToFlash, "source")):
    if ".svn" in dirs:
      dirs.remove(".svn")
    for f in files:
      if f.upper().endswith(".F90"):
        f90Files.append(os.path.normpath(os.path.join(root, f)))

  if not f90Files:
    print "No F90 files found."
    sys.exit(1)

  interfacesDict = InterfacesDict(pathToFlash)

  unitnames = []
  if "-auto" in args:
    unitnames.extend(interfacesDict.getAllInterfacedUnits())
    args.remove("-auto")

  if len(args) > 1:
    unitnames.extend(args[1:])

  adjust(f90Files, interfacesDict, unitnames, dryRun)
