#!/usr/bin/env python
import sys, os, re

# matches a robodoc 'function' or 'header' header,
# with or without an 'internal' marker
fHeaderPat = re.compile("(\!\!)\*\*\*\*i?([f|h])\*")
cHeaderPat = re.compile("(/)\*\*\*\*i?([f|h])\*")

acceptedExtensions  = [(".F90", fHeaderPat), (".c", cHeaderPat)]
excludedDirectories = [".svn"]

def hasAcceptedExtension(f):
  """
  if list comprehension returns list of len > 0, it means
  'f' ends in one of the accepted extensions
  """
  for acceptedExtension in acceptedExtensions:
    if f.endswith(acceptedExtension[0]):
      return acceptedExtension[1]  # return the regex that matches a robodoc
                                   # header for a file with this extension
  return None

def isExcludedDirectory(d):
  """
  if list comprehension returns list of len > 0, it means
  'd' is one of the excluded directories
  """
  return len([1 for dir in excludedDirectories if d == dir])

def fixHeaders():
  """
  Recursively search through all files of appropriate extensions in cwd.
  If a robodoc header is found, determine the correct path and filename
  that should follow that header. If this is not the path and filename
  actually present, correct it.

  A recognized robodoc header (on a Fortran file) will have the form:

        !!****f*    or    !!****h*
  """
  cwd = os.getcwd()
  items = os.listdir(".")

  for item in items:
    if os.path.isdir(item) and not isExcludedDirectory(item):
      os.chdir(item)
      fixHeaders()
      os.chdir(cwd)
    elif os.path.isfile(item):
      headerPat = hasAcceptedExtension(item)
      if headerPat:
        firstLine = open(item).readline().strip()
        m = headerPat.match(firstLine)
        if m:
          # 'shortPathToDir' will be used to help determine 'targetLine' below.
          # The path and filename that we want written into the robodoc header
          # should assume an uppermost root directory at the top of a working
          # copy of the FLASH code, not the root of the local machine as passed
          # in from 'cwd'. Therefore we truncate 'pathToFlash'
          shortPathToDir = cwd[len(pathToFlash):].strip("/")

          # Since different languages have different ways of notating a comment,
          # robodoc headers also look different in the different languages. The
          # first group of the regex captures this comment mark for the different
          # languages.
          commentMark = m.group(1)

          # To determine whether this function should be marked "internal", we
          # find the first instance in 'shortPathToDir' of a slash followed by a
          # capital letter. If there are more slashes beyond this point, then the
          # function is not immediately under the first directory starting with a
          # capital, and must be marked with an "i".
          m2 = re.search("[/][A-Z]", shortPathToDir)
          if m2 and shortPathToDir.find("/", m2.end()) < 0:
            internalInfix = ""
          else:
            internalInfix = "i"

          # 'targetLine' is what we *should* find as the first line of this file:
          # the robodoc header, one space, and the filename minus the extension
          # NB! Truncating the extension by using rfind will only work if the
          # filename has a "." in it somewhere, but we know it must because it
          # passed the 'hasAcceptedExtension' test earlier.
          targetLine = (commentMark + "****" + internalInfix + m.group(2) + "* " +
                        os.path.join(shortPathToDir, item[0:item.rfind(".")]))
          if targetLine != firstLine:
            # If the file's first line is wrong, replace it with 'targetLine',
            # and keep the rest the same
            text = open(item).read()
            open(item, "w").write(targetLine + text[text.find("\n"):])
            #print "corrected header for: " + os.path.join(cwd, item)

# this script assumes it is being run from inside FLASH3/tools/scripts
cwd = os.getcwd()
pathToFlash = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0]))))
os.chdir(os.path.join(pathToFlash, "source"))
fixHeaders()
os.chdir(cwd)
