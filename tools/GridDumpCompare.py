#!/usr/bin/env python
"""
This is the python version of Anshu's GridDumpCompare.F90. The utility
compares two files written by the function Grid_dump in FLASH and
determines if they are identical.

The program is invoked from the command line with paths to the two files
to be compared as arguments:

  ./GridDumpCompare.py <pathToFile1> <pathToFile2>
"""
import array, os, sys

def main(twoFiles, tol):
  arraysFromFiles = []

  for i in range(2):
    if not os.path.isfile(twoFiles[i]):
      print "ERROR: \"%s\" does not exist or is not a file." % twoFiles[i]
      sys.exit(1)
    else:
      newArray = array.array("d")  # make array for "d"ouble size reals
      numItems = os.stat(twoFiles[i])[6]/newArray.itemsize
      newArray.fromfile(open(twoFiles[i], "rb"), numItems)
      arraysFromFiles.append(newArray)

  if len(arraysFromFiles[0]) != len(arraysFromFiles[1]):
    print "ERROR: the arrays in the two files are different sizes, and cannot be compared."
    print "The array from file \"%s\" contains %s items." % (twoFiles[0], len(arraysFromFiles[0]))
    print "and the array from file \"%s\" contains %s items." % (twoFiles[1], len(arraysFromFiles[1]))
    sys.exit(1)

  # else
  maxDiff = 0.0
  diffLoc = -1
  for i in range(len(arraysFromFiles[0])):
    if abs(arraysFromFiles[0][i] - arraysFromFiles[1][i]) >= maxDiff:
      maxDiff = abs(arraysFromFiles[0][i] - arraysFromFiles[1][i])
      diffLoc = i

  if maxDiff == 0.0:
    print "The two files are identical."
    sys.exit(0)
  elif maxDiff < tol:
    print "The two files are similar within the given tolerance."
    sys.exit(0)
  else:
    print "The two files show a maximum difference of %s" % maxDiff
    print "at word number %s" % diffLoc
    print "where \"%s\" has a value of %s" % (twoFiles[0], arraysFromFiles[0][diffLoc])
    print "and   \"%s\" has a value of %s" % (twoFiles[1], arraysFromFiles[1][diffLoc])
    sys.exit(1)


def usage():
  print "usage:"
  print "  GridDumpCompare.py [-t tolerance] <pathToFileA> <pathToFileB>"
  sys.exit(0)

if __name__=="__main__":
  args = sys.argv
  if args.count("-t") > 0:
    tolIndex = args.index("-t")
    if len(args) == tolIndex + 1:
      usage()
    # else
    try:
      tol = float(args[tolIndex + 1])
    except:
      usage()
    # else
    del args[tolIndex:tolIndex+2] # this gets rid of exactly two arguments
  else:
    tol = 0.0

  if len(sys.argv) != 3:
    usage()
  else:
    main(args[1:], tol)
