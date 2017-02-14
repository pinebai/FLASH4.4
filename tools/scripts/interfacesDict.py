import os, re

__all__ = ["INTERFACE_FILE_MISSING", "INTERFACE_DECL_IN_ORDER", "INTERFACE_DECL_NO_FILE",
           "INTERFACE_DECL_MISSING", "INTERFACE_DECL_EXTRA", "INTERFACE_CASE_MISMATCH",
           "InterfacesDict"]

INTERFACE_FILE_MISSING  = 10
INTERFACE_DECL_IN_ORDER = 20
INTERFACE_DECL_NO_FILE  = 30
INTERFACE_DECL_MISSING  = 40
INTERFACE_DECL_EXTRA    = 50
INTERFACE_CASE_MISMATCH = 60

class InterfacesDict(dict):
  """
  This class represents a dictionary whose keys are paths to API-level
  directories (e.g. source/monitors/Timers) and whose values are either:
    
    INTERFACE_FILE_MISSING if no "<unit-name>_interface.F90" file
    at that location
    
  or:
    
    another dictionary whose keys are the names of the subroutines found
    in the interface file and whose values are one of:
      INTERFACE_DECL_IN_ORDER - all is well
      INTERFACE_DECL_NO_FILE  - this subroutine needs no API-level file
      INTERFACE_DECL_MISSING  - file at API-level has no interface
      INTERFACE_DECL_EXTRA    - interface matches no API-level file
      INTERFACE_CASE_MISMATCH - interface differs from API-level
                                file only in case of its name
  """
  doubleQuotesPat = re.compile(r'(?<!\\)".*?(?<!\\)"', re.IGNORECASE)
  singleQuotesPat = re.compile(r"(?<!\\)'.*?(?<!\\)'", re.IGNORECASE)
  commentPat      = re.compile(r"^!.*$", re.IGNORECASE)
  #interfacePat   = re.compile(r"^\s*interface\s+(?P<intname>[-\w]+)\s*$", re.IGNORECASE)
  interfacePat    = re.compile(r"^\s*interface(?:\s+(?P<intname>[-\w]+))?\s*$", re.IGNORECASE)
  endInterfacePat = re.compile(r"^\s*end\s+interface\s*$")
  subroutinePat   = re.compile(r"^\s*(?:recursive\s+)?subroutine\s*(?P<subname>[-\w]+)\s*\(", re.IGNORECASE)

  def stripQuotesAndComments(self, line):
    line = self.doubleQuotesPat.sub("", line)
    line = self.singleQuotesPat.sub("", line)
    line = self.commentPat.sub("", line)
    return line

  def __init__(self, pathToFlash):
    """
    constructor
    """
    self.unitToInterfaceMap = {}
    self.errorsFound = False

    cwd = os.getcwd()
    os.chdir(pathToFlash)

    pathsToAPILevelDirs = []
    # walk source tree adding the first capitalized
    # directory to 'pathsToAPILevelDirs'
    for root, dirs, files in os.walk("source"):
      i=0
      dirs.sort()
      while i < len(dirs):
        if dirs[i][0].isupper():
          pathsToAPILevelDirs.append(os.path.join(root, dirs[i]))
          del dirs[i]
        else:
          i+=1

    for pathToAPILevelDir in pathsToAPILevelDirs:
      APILevelUnit = os.path.basename(pathToAPILevelDir)
      interfaceFilename = APILevelUnit + "_interface.F90"
      pathToInterfaceFile = os.path.join(pathToAPILevelDir, interfaceFilename)

      # Make a note if expected interface file isn't there
      if not os.path.isfile(pathToInterfaceFile):
        self[pathToAPILevelDir] = INTERFACE_FILE_MISSING
        self.errorsFound = True
        continue

      # else map this unit's name to its location in "source"
      self.unitToInterfaceMap[APILevelUnit] = pathToAPILevelDir

      # Parse the file to see which interfaces it lists
      self[pathToAPILevelDir] = {}
      subnamesDict = self[pathToAPILevelDir]  # get a handle that's easier to reference

      lines = open(pathToInterfaceFile).read().split("\n")
      i = 0
      while i < len(lines):
        # strip out all double quotes, single quotes, and comments
        lineCopy = self.stripQuotesAndComments(lines[i])
        if lineCopy:
          m0 = self.interfacePat.match(lineCopy)
          if m0:
            intname = m0.group("intname")
            while not self.endInterfacePat.match(lineCopy):
              i+=1
              lineCopy = self.stripQuotesAndComments(lines[i])
              m1 = self.subroutinePat.match(lineCopy)
              if m1:
                subname = m1.group("subname")
                subnamesDict[subname] = INTERFACE_DECL_IN_ORDER  # may change in a couple of lines
                if intname and subname != intname:  # There may be no 'intname' if the interface
                                                    # only wraps a single subroutine

                  # This subroutine is inside an interface whose name differs from
                  # its own; e.g. "Timers_stopIndex" is under "Timers_stop". Such a
                  # subroutine needs no corresponding stub file at the API level.
                  subnamesDict[intname] = INTERFACE_DECL_IN_ORDER
                  subnamesDict[subname] = INTERFACE_DECL_NO_FILE
        i+=1

      # Note all files in 'pathToAPILevelDir' that start with 'APILevelUnit'
      # (unless followed by "_interface" or "_data")
      APIFilePat       = re.compile(r"(%s(?:_(?!interface|data)\w+)?)\.F90$" % APILevelUnit)
      APIStubFilenames = []

      # Get list of all files under 'pathToAPILevelDir'
      allFilesUnderPath = [item for item in os.listdir(pathToAPILevelDir) if
                           os.path.isfile(os.path.join(pathToAPILevelDir, item))]

      for f in allFilesUnderPath:
        m = APIFilePat.match(f)
        if m:
          APIStubFilenames.append(m.group(1))

      for APIStubFilename in APIStubFilenames:
        if APIStubFilename not in subnamesDict.keys():
          correctlyCasedSubname = self.getCorrectlyCasedSubname(APIStubFilename)
          if correctlyCasedSubname:
            # There is a subroutine in this Unit that corresponds to
            # 'APIStubFilename', but with a different capitalization pattern.
            subnamesDict[correctlyCasedSubname] = INTERFACE_CASE_MISMATCH
          else:
            subnamesDict[APIStubFilename] = INTERFACE_DECL_MISSING
          self.errorsFound = True

      for declaredInterface in subnamesDict.keys():
        # The interface could actually be a subroutine inside an interface
        # like "Timers_stopIndex" under "Timers_stop". If that's the case,
        # we don't need to look for a corresponding stubfile.
        if subnamesDict[declaredInterface] == INTERFACE_DECL_NO_FILE:
          continue
        # else
        if declaredInterface not in APIStubFilenames:
          # It won't be in here if it's already been declared 'case-mismatch'
          if subnamesDict[declaredInterface] == INTERFACE_CASE_MISMATCH:
            continue
          else:
            subnamesDict[declaredInterface] = INTERFACE_DECL_EXTRA
            self.errorsFound = True

    os.chdir(cwd)


  def getPathAndSubroutinesForUnit(self, unitname):
    """
    Return a tuple consisting of a path from "source" to the directory
    containing the _interface.F90 file that corresponds to 'unitname',
    and the dictionary representing the contents of that file.
    """
    unitname = self.getCorrectlyCasedUnitname(unitname)
    if self.unitToInterfaceMap.has_key(unitname):
      pathToAPILevelDir = self.unitToInterfaceMap[unitname]
      return (pathToAPILevelDir, self[pathToAPILevelDir])
    return (None, None)

  def getAllInterfacedUnits(self):
    """
    Return list of all API-level unitnames for which an _interface.F90
    file exists. These can be used to retrieve other information from
    getPathAndInterfaceForUnit().
    """
    return self.unitToInterfaceMap.keys()

  def getCorrectlyCasedUnitname(self, unitname):
    """
    The submitted unit name ('unitname') will be compared in a
    case-insensitive comparison with all unitnames known to this
    instance of InterfaceDict. If a match is found, return the
    unit name with correct capitalization, where correct is defined
    as whatever is found in the _interface.F90 file.
    """
    allUnitnames = self.getAllInterfacedUnits()
    upperToCasedDict = dict([(casedUnitname.upper(), casedUnitname) for casedUnitname in allUnitnames])
    return upperToCasedDict.get(unitname.upper(), None)

  def getCorrectlyCasedSubname(self, subname):
    """
    The submitted subroutine name ('subname') will be compared in a
    case-insensitive comparison with all subroutine names valid for
    its unit, which is derived from its (presumably extant) prefix.
    If a match is found, return the subroutine name with correct
    capitilazation, where correct is defined as whatever is found in
    the _interface.F90 file.
    """
    unitname = subname.split("_",1)[0]
    pathToAPILevelDir, subnamesDict = self.getPathAndSubroutinesForUnit(unitname)
    if not pathToAPILevelDir:
      return None
    # else
    upperToCasedDict = dict([(casedSubname.upper(), casedSubname) for casedSubname in subnamesDict])
    return upperToCasedDict.get(subname.upper(), None)
