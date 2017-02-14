#!/usr/bin/env python
import sys, os, re

# regexps used to pare down a line of text to essential elements
doubleQuotesPat = re.compile(r'(?<!\\)".*?(?<!\\)"')
singleQuotesPat = re.compile(r"(?<!\\)'.*?(?<!\\)'")
commentPat      = re.compile("!.*$")

class SubroutineInfo:
  def __init__(self):
    self.type                  = ""  # subroutine or function
    self.subname               = ""  # name of this subroutine/function
    self.subroutineLineIndex   = -1  # line number of this subroutine's "subroutine" declaration
    self.poundIncludeLineIndex = -1  # line number of a '#include' that should come out
    self.useOnlyLine           = ""  # text of a pat2 'use', "&"-segments rejoined if necessary
    self.useOnlyLineIndices    = []  # line numbers over which a pat2 'use' appears
    self.indent                = -1  # width of whitespace before pat2 'use'
    self.otherUseLineIndex     = -1  # line number where the first non-pat2 'use' appears
    self.useOnlySubnames       = []  # subroutines already listed in the pat2 'use'
    self.calledSubnames        = []  # list of calls to subroutines in body of code that match 'pat8'.

def stripQuotes(line):
  line = doubleQuotesPat.sub("", line)
  line = singleQuotesPat.sub("", line)
  return line

def stripComments(line):
  line = commentPat.sub("", line)
  return line

def stripQuotesAndComments(line):
  line = stripQuotes(line)
  line = stripComments(line)
  return line

def getUseAndCallInfo(unitName, lines):
  """
  Return various data related to unit 'unitName' as found in the
  file whose lines are captured in the list 'lines'. This module
  is useable as a standalone or via "adjustInterfaceCall.py"
  """
  pat0 = re.compile("^\s*(?:(?:recursive|real|integer|logical)\s+)?(program|subroutine|function)\s*([-\w]*)\s*", re.IGNORECASE)
  pat1 = re.compile("^\s*#\s*include\s*\"(%s_interface\.h)\"" % unitName, re.IGNORECASE)
  pat2 = re.compile("^(\s*)(use\s+?%s_interface\s*?(?:,\s*?ONLY\s*:)?)" % unitName, re.IGNORECASE)
  pat3 = re.compile("^(\s*)(use .*)", re.IGNORECASE)
  pat4 = re.compile("^\s*interface(?:\s|$)", re.IGNORECASE)
  pat5 = re.compile("^\s*end\s+interface", re.IGNORECASE)
  pat6 = re.compile("^\s*end\s+(program|subroutine|function)", re.IGNORECASE)
  pat7 = re.compile("%s(?:_\w+)?(?:\s|,|$)" % unitName, re.IGNORECASE)  # *only* for finding fns listed after the pat2 'use'
  pat8 = re.compile("(?:^|.*?\s*)(call\s+(%s(?:_\w+)?))(\w+)?\s*" % unitName, re.IGNORECASE)

  allSIs = []
  si = SubroutineInfo()

  i = 0
  while i < len(lines):
    lineCopy = lines[i]
    if not pat1.search(lineCopy):  # pound-includes are followed by the unit-name in quotes
      lineCopy = stripQuotes(lineCopy)  # strip out strings in single or double quotes
    lineCopy = stripComments(lineCopy)  # strip everything to right of comment
    lineCopy = lineCopy.rstrip()

    if lineCopy:
      m0 = pat0.search(lineCopy)  # find type/name of this subroutine/function
      m1 = pat1.search(lineCopy)  # find a "#include <unitName>_interface.h" that should come out
      m2 = pat2.search(lineCopy)  # find already-extant instance of "use" for our unit name
      m3 = pat3.search(lineCopy)  # find already-extant instance of "use" for something else
      m4 = pat4.search(lineCopy)  # find whether we've entered an "interface" block
      m6 = pat6.search(lineCopy)  # find end of a subroutine block

      if m0:
        # this is a 'subroutine' or 'function' declaration. The name should follow
        si.type, si.subname = m0.groups()
        while lines[i].strip().endswith("&"):
          i+=1
        si.subroutineLineIndex = i+1  # actually the first line *under* the subroutine
                                      # declaration where it's safe to insert something
      elif m1:
        # this line contains a pound-include of the kind we're trying to eliminate
        # PURE HACK to exclude "Eos.h" and others which are not actually interface
        # files the way "Logfile.h" is.
        # KW 2009-06-30: This hacky if is not necessary any more, since now we are only trying to
        # eliminate includes of the form "#include <unitName>_interface.h". Leaving it in place
        # anyway, it cannot hurt.
        if m1.group(1) != "Eos.h" and m1.group(1) != "MHD.h" and m1.group(1) != "Multispecies.h":
          si.poundIncludeLineIndex = i

      elif m2:
        if len(si.useOnlyLineIndices) == 0:
          # collect all pieces of this line, which may be broken
          # over several lines via the "&" character.
          si.indent = len(m2.group(1))
          while lines[i].strip().endswith("&"):
            si.useOnlyLine += (lines[i].strip().rstrip("&") + " ")
            si.useOnlyLineIndices.append(i)
            i+=1
          si.useOnlyLine += lines[i].strip()
          si.useOnlyLineIndices.append(i)

          # chop off the part of the string that actually
          # has the 'use interface, ONLY:' part
          si.useOnlyLine = si.useOnlyLine[len(m2.group(2)):]
          si.useOnlySubnames = [item.rstrip(",") for item in pat7.findall(si.useOnlyLine)]
        else:
          while lines[i].strip().endswith("&"):
            i+=1

      elif m3:
        if si.otherUseLineIndex == -1:
          # This line is a 'use' line, but it doesn't match pat2. We
          # might need this if we need to insert a brand-new pat2 'use'
          # (as opposed to replacing an already existing one)
          if len(si.useOnlyLineIndices) == 0:  # only if we haven't already set this in the 'm2' case
            si.indent = len(m3.group(1))
          while lines[i].strip().endswith("&"):
            i+=1
          si.otherUseLineIndex = i+1  # actually the first line *under* the first 'use'
                                      # instance where it's safe to insert something
        else:
          while lines[i].strip().endswith("&"):
            i+=1

      elif m4:
        # We've entered an "interface" declaration, which can contain
        # its own "subroutine" declarations and confuse the script.
        # So keep skipping lines until we're out of it.
        while not pat5.match(lines[i]):
          if i < len(lines):
            i += 1

      elif m6:
        if m6.group(1).upper() == si.type.upper():
          # We've reached the end of a subroutine block, so close
          # out this instance of SubroutineInfo and start a new one.
          allSIs.append(si)
          si = SubroutineInfo()

      else:
        # only take a value for indent if we're past the 'subroutine'
        # declaration and no 'indent' value has been set before now.
        if (si.subname and len(lines[i].strip()) > 0 and si.indent == -1):
          si.indent = len(re.match("^(\s*)", lines[i]).group())

        # Capture all calls to API-level subroutines in 'si.calledSubnames'
        # Each element in the list is a tuple containing the subroutine name,
        # the complete call to the subroutine, and the line on which the call
        # occurs. 'pat8' captures the entire "call" context so that later when
        # we check capitalization consistency between these calls and the names
        # of subroutines at the API-level, false "corrections" will not be made
        # to other instances of the subroutine name that occur outside the call
        # statement.
        while lines[i].strip().endswith("&"):
          lineCopy = stripQuotesAndComments(lines[i]).strip()
          si.calledSubnames.extend([(m[1], m[0], i) for m in pat8.findall(lineCopy) if not m[2]])
          i+=1
        lineCopy = stripQuotesAndComments(lines[i]).strip()
        si.calledSubnames.extend([(m[1], m[0], i) for m in pat8.findall(lineCopy) if not m[2]])

    i+=1

  return allSIs


if __name__=="__main__":
  if len(sys.argv) < 2:
    print "You must supply a unit-name"
    sys.exit(1)
  elif len(sys.argv) < 3:
    print "You must supply a unit-name and a filename"
    sys.exit(1)
  else:
    unitName = sys.argv[1]
    lines = open(sys.argv[2]).read().split("\n")
    allSIs = getUseAndCallInfo(unitName, lines)
    for si in allSIs:
      print "subroutineName is: %s" % si.subname
      print "subroutineLineIndex is: %s" % si.subroutineLineIndex
      print "poundIncludeLineIndex is: %s" % si.poundIncludeLineIndex
      print "useOnlyLineIndices is: %s" % si.useOnlyLineIndices
      print "otherUseLineIndex is: %s" % si.otherUseLineIndex
      print "indent is: %s" % si.indent
      print "useOnlySubnames is: %s" % si.useOnlySubnames
      print "calledSubnames is: %s" % si.calledSubnames
      print ""
