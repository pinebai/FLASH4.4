
# Code which parses the command line arguments 

# all functions in parseCmd must be explicitly referred to
__all__ = []

########################## Internal implementation ##########################

import string, re, getopt, sys, os
import globals, lazyFile
from globals import *
from utils import *

# default arguments
ADDL_DEF_ARGS = ["+default"]

# options without arguments
WITHOUT_ARGS = ["auto","1d","2d","3d","portable",
                "makehide", "curvilinear",
                "opt","debug","test", "index-reorder", "strictparams",
                "fbs","help", "noclobber", "nofbs"] # nofbs disables fixedBlockSize
# options with arguments
WITH_ARGS = ["maxblocks","nxb","nyb","nzb","verbose","site","ostype",
             "defines", "objdir","with-unit", "unit", "with-library",
             "without-unit","without-library","kill-unit","unitsfile",
             "makefile", "library", "datafiles", "parfile", "tau",
             "gridinterpolation", "geometry", "particlemethods"]


USAGE="""usage:  setup <problem-name> [options] [VAR=VALUE]...

   problem-name: see source/Simulation/SimulationMain directory
   options: 

   (Science Options)
            -auto -[123]d 
            -maxblocks=<#> -nxb=<#> -nyb=<#> -nzb=<#>
            -with-unit=<unit> -with-library=<libname>[,args]
            -without-unit=<unit> -without-library=<libname>

   (Setup and Make Options)
            -verbose=[DEBUG|INFO|WARN|IMPINFO|ERROR] 
            [-site=<site> | -ostype=<ostype>] 
            -makefile=<extension>
            [-opt| -debug | -test ] 
            -objdir=<relative obj directory> 
            -defines=<defines> -unitsfile=<filename>
            -datafiles=<wildcard> -parfile=<filename>
            -fbs -nofbs -tau=<makefile>

   (Misc Options)
            -makehide -noclobber -portable -help

   * For GNU compatibility, options may be prefixed by -- instead of - as well
   * -unit and -library are considered equivalent to 
     -with-unit and -with-library respectively.
   * For information regarding the [VAR=VALUE] options and using 'setup variables' 
     refer to User's Guide.  
   * To read how shortcuts work see README.shortcuts in your bin directory
"""

def usage():
    """Print usage info and exit"""
    print USAGE
    ans = raw_input("\nDo you want to see a list of shortcuts I know about [Y/n]?")
    ans = ans.replace("\n","")
    if ans.lower() in ["y","yes",""]:
       sitems = GVars.shortcuts.items()
       sitems.sort()
       c = 0
       for (k,v) in sitems:
           c = max(c,len(k))
       tpl = "   %%-%ds %%s" % (c+2) # template for printing the dictionary
       print "\nTo use a shortcut add '+shortcut' to your setup line."
       print "For example ./setup Sod -auto +ug\n"
       for k,v in sitems:
           print tpl % (k," ".join(v))
    raise SetupError("")

def cleanupCommandLine(args):
   """Return a lists where all single - options are converted to --"""
   ans = []
   for m in args:
       if len(m) > 1 and m[0] == "-" and m[1] != "-":
             ans.append("-"+m)
       else: ans.append(m)
   return ans

# processes an existing file and returns a dictionary of shortcuts
def getShortcutDict(filename):
   sfd = file(filename)
   shortcuts = {}
   GVars.out.put("Processing Shortcut file: %s" % filename,globals.IMPINFO)
   for line in sfd:
      if not line: continue
      if line[0] == "#": continue
      line = line.strip()
      if line and line[-1] in ["\r","\n"]: line = line[:-1]
      if line and line[-1] in ["\r","\n"]: line = line[:-1]
      if not line: continue
      # now we have a valid line
      parts = line.split(":")
      value = cleanupCommandLine([x.strip() for x in parts[1:] if x])
      shortcuts[parts[0].lower()] = value
   return shortcuts

# shortcuts is a dictionary mapping shortcut to list of arguments
# shortcut X is invoked as "+x" option to setup
def getShortcuts():
   if not os.environ.has_key("SETUP_SHORTCUTS"):
      GVars.out.put("No Shortcut file specified using default",globals.DEBUG)
      sfiles = ["setup_shortcuts.txt"]
   else:
      sfiles = os.environ["SETUP_SHORTCUTS"].split(":")
   shortcuts = {}
   for name in sfiles:
      sfile = os.path.abspath(name)
      if not os.path.exists(sfile):
         GVars.out.put("Unable to open %s. Ignoring" % sfile,globals.WARN)
      else: shortcuts.update(getShortcutDict(sfile))
   return shortcuts

def expandShortcuts(args,RecLimit=255):
   """Expand all shortcuts specified in args. 
      No more than RecLimit expansions should be required. 
      This limit can be used to detect circular references"""
   ans = []
   count = RecLimit
   GVars.shortcuts = getShortcuts()
   shortcuts = GVars.shortcuts
   while args:
      scut = args[0]
      del args[0]
      if scut[0] != globals.SHORTCUT_CHAR: 
         ans.append(scut) 
      else:
         # found a shortcut
         cand = []
         # find a shortcut starting with given letters (case insensitive)
         for x in shortcuts.keys():
             y = x.lower()
             if y.startswith(scut[1:].lower()): cand.append(x)
         if len(cand) == 1:
            scut = cand[0] # only one candidate
         elif shortcuts.has_key(scut[1:].lower()): # specified key as such exists
            scut = scut[1:].lower() # pick given key
         # now process the shortcut
         if not shortcuts.has_key(scut): # invalid shortcut
            GVars.out.put("\n***WARNING*** Ignoring unknown shortcut %s while expanding %s.\n" % (scut," ".join(args)),globals.IMPINFO)
         elif count > 0:
            count = count - 1
            args[0:0] = shortcuts[scut] # insert expansion of shortcut in front of args
         else:
            raise SetupError("Too many shortcuts. Circular reference?")
   return ans

def setSetupVars(strlist):
   """Given a list of strings. Identifies those of the form A=B and 
     updates GVars.setupVars.  Returns unprocessed list of strings"""
   ans = []
   # declare setup variables for all top level units
   # these are boolean variables which decide whether
   # to include this unit or not
   for x in GVars.topUnitNames: # by default do not include any units
       GVars.setupVars.addunit(x,False)
   for opt in strlist:
       p = opt.find("=")
       if p < 0: # no equal found
          ans.append(opt)
       else:
          GVars.setupVars.set(opt[:p],opt[p+1:]) # set variable to value
   return ans
 
def custom_getopt(longoptions):
   """Parse command line arguments, working in gnu mode using only longoptions.
     GNU MODE =  allow non-option arguments to be followed by options.
     We also allow user to use shortcuts as well as define variables like make"""
   max_shortcuts = 255 # no more than 255 shortcuts allowed per invocation
   optvallist = []
   rest = []
   toparse = cleanupCommandLine(ADDL_DEF_ARGS+sys.argv[1:]) # double "-" to "--"
   fulllist = expandShortcuts(toparse) # expand the shortcuts
   fullCmdList = [x for x in fulllist if x and (x[0] == "-" or x.find("=") < 1) ] # remove non-options with "=" in it
   GVars.fullCmdLine = " ".join(fullCmdList) # this gets printed
   toparse = fulllist[:] # make a copy of the list
   while toparse:
    optval_sublist, toparse = getopt.getopt(toparse,"",longoptions)
    optvallist.extend(optval_sublist)
    if toparse: 
       badoption = toparse[0]
       del toparse[0]
       rest.append(badoption)
   return (fulllist,optvallist,setSetupVars(rest))

def parseCommandLine():
    # process all the options
    # those requiring arguments suffixed with =
    longopts = WITHOUT_ARGS + [x+"=" for x in WITH_ARGS]

    try:
       (fullcmdline,optvallist,rest) = custom_getopt(longopts)
    except getopt.GetoptError,e:
       GVars.out.put(str(e),globals.ERROR)
       usage() # print usage info

    # Not given any program name
    if not rest: usage()

    # if given many the last one wins
    if len(rest) > 1:
       GVars.out.put("         WARNING: Multiple problem names given. winner="+rest[-1],globals.WARN)
    GVars.simulationName = rest[-1] # store program name

    # all acceptable values for "--verbose" keyword
    vrblevels = {"DEBUG":globals.DEBUG,
                 "PPDEBUG":globals.PPDEBUG,
                 "PPWARN":globals.PPWARN,
                 "WARN":globals.WARN,
                 "INFO":globals.INFO,
                 "IMPINFO":globals.IMPINFO,
                 "ERROR":globals.ERROR}

    # all acceptable values for "--gridinterpolation" keyword
    allGridInterpolations = {"MONOTONIC":globals.GRID_INTERP_MONOTONIC,
                             "NATIVE":globals.GRID_INTERP_NATIVE}

    # all acceptable values for "--geometry" keyword
    allGeometries = {"CARTESIAN":globals.GRID_GEOM_CARTESIAN,
                     "CYLINDRICAL":globals.GRID_GEOM_CYLINDRICAL,
                     "SPHERICAL":globals.GRID_GEOM_SPHERICAL,
                     "POLAR":globals.GRID_GEOM_POLAR}

    withUnits = {} # dictionary to store list of units to add (prevents duplicate entries)

    for (arg,val) in optvallist:
        if   arg == '--portable':         GVars.portable = 1
        elif arg == '--noclobber':        GVars.noClobber = 1
        elif arg == '--auto':             GVars.auto = 1
        elif arg == '--datafiles':        GVars.datafiles.append(val)
        elif arg == '--parfile':          GVars.parfile = val
        elif arg == '--debug':            GVars.buildFlag = "DEBUG"
        elif arg == '--test':             GVars.buildFlag = "TEST"
        elif arg == '--opt':              GVars.buildFlag = "OPT"
        elif arg == '--maxblocks':          GVars.maxblocks = int(val)
        elif arg == '--makehide':           GVars.makedisplay = 0
        elif arg == '--makefile':           GVars.makefileext = "."+val
        elif arg == '--nxb':                GVars.nxb = int(val)
        elif arg == '--nyb':                GVars.nyb = int(val)
        elif arg == '--nzb':                GVars.nzb = int(val)
        elif arg == '--site':               GVars.build_site = val
        elif arg == '--tau':                GVars.build_tau = val
        elif arg == '--ostype':             GVars.build_os = val
        elif arg == '--objdir':             GVars.objectDir = val
        elif arg == '--help':               usage()
        elif arg == '--unitsfile':          GVars.unitsFile = val
        elif arg == '--fbs':                GVars.setupVars.set("fixedBlockSize",True)
        elif arg == '--nofbs':              GVars.setupVars.set("fixedBlockSize", False)
        elif arg == '--strictparams':       GVars.strictParams = 1
        # DEV 'curvilinear'
        # * originally used to have the same effect as -gridinterpolation=monotonic
        #   does now, in addition to #defining GRID_CURVILINEAR 1 in Flash.h.
        # * Then there came a time when
        #   o  -curvilinear was deprecated, and
        #   o  -gridinterpolation=monotonic had the effect of always #defining
        #      GRID_CURVILINEAR 1 in Flash.h.
        # * Now the time has come when
        #   o  -gridinterpolation=monotonic does NOT automatically #define GRID_CURVILINEAR;
        #   o  (but -geometry={cylindrical,spherical,polar} still does;)
        #   o  -geometry={cylindrical,spherical,polar} implies -gridinterpolation=monotonic;
        #   o  -curvilinear is not deprecated any more, since it can be used to force
        #      #defining GRID_CURVILINEAR 1 in Flash.h when this is not automatically
        #      done any more. That is, -curvilinear can be useful (and has only an effect)
        #      when -geometry={cylindrical,spherical,polar} is not requested.
        #      However, this is rarely useful..
        elif arg == '--curvilinear':
             GVars.out.put("\n***************************** INFO *********************************\n" +
                           "The -curvilinear flag is nearly always unnecessary.\n" +
                           "Use --geometry={cylindrical,spherical,polar} to configure\n" +
                           "FLASH for a specific non-Cartesian geometry.\n" +
                           "(These flags are not needed at all when using PARAMESH in LIBRARY\n" +
                           "MODE, i.e., ParameshLibraryMode=True, or when using Paramesh4dev.)\n" +
                           "********************************************************************\n",globals.WARN)
             GVars.curvilinear = 1
        elif arg in ["--1d","--2d","--3d"]: 
             GVars.dimension = int(arg[2]) # i.e. 1 or 2 or 3
             GVars.setupVars.set("nDim",GVars.dimension)
        elif arg == "--index-reorder":            
             GVars.setupVars.set("GridIndexReordered",True) # for use in Config file
             GVars.indexReorder = True # for use in Makefile
             GVars.defines.append("-DINDEXREORDER") # for use in code
             GVars.definesNames.append("-DINDEXREORDER")
        elif arg in ["--with-unit","--unit"]:
             if val.endswith(os.sep): val = val[:-1]
             if val.startswith("source"+os.sep): val = val[7:]
             withUnits[val] = 1
             GVars.setupVars.addunit(val,True) # set corresponding setup var to true
        elif arg == "--without-unit": # remove all units added which come under specified unit
             # also add to list of units to be ignored when handling REQUESTS keyword
             if val.endswith(os.sep): val = val[:-1] 
             rmlist = [x for x in withUnits.keys() if x == val or x.startswith(val+os.sep)]
             for x in rmlist:
                 del withUnits[x]
             GVars.withoutUnits[val] = 1
        elif arg == "--kill-unit":
             GVars.killUnits[val] = 1
        elif arg in ["--with-library","--library"]: 
             if not val: continue # no name given --> ignore 
             parts = val.split(",")
             libname = parts[0]
             args = string.join(parts[1:]," ") # replace commas with space
             if len(args) >= 2 and args[0] == args[-1] and args[0] in ['"',"'"]: # argument has been quoted
                args = args[1:-1]
             GVars.withLibraries[libname.lower()] = args
        elif arg == "--without-library":
             val = val.lower()
             if GVars.withLibraries.has_key(val): del GVars.withLibraries[val]
             GVars.withoutLibraries[val] = 1
        elif arg == '--verbose': # set verbosity level
             if not val: continue # no argument dont change level
             if vrblevels.has_key(val.upper()): 
                GVars.verbose = vrblevels[val.upper()]
             else:
                GVars.out.put("Unrecognized verbosity level [%s]" % val,globals.ERROR)
                usage()
        elif arg == '--gridinterpolation':  # set grid interpolation
             if not val: continue  # no argument; don't change interpolation
             if allGridInterpolations.has_key(val.upper()):
                 GVars.gridInterpolation = allGridInterpolations[val.upper()]
             else:
                 GVars.out.put("Unrecognized grid interpolation [%s]" % val, globals.ERROR)
                 usage()                 
        elif arg == '--geometry': # set geometry
             if not val: continue # no argument; don't change geom
             if allGeometries.has_key(val.upper()):
                 GVars.gridGeometry = allGeometries[val.upper()]
                 if val.upper() != "CARTESIAN":
                     # All geometries other than cartesian (i.e. cylindrical,
                     # spherical, and polar) are automatically curvilinear.
                     # (See note above on -curvilinear.)
                     GVars.curvilinear = 1
             else:
                 GVars.out.put("Unrecognized geometry [%s]" % val, globals.ERROR)
                 usage()
        elif arg == "--defines": # declare additional CPP/FPP stuff
             if not val: # kill existing defines
                GVars.defines=[]
                GVars.definesNames=[]
                continue 
             for x in val.split(","):
                 x = x.strip()
                 p = x.find("=")
                 if p < 0: # no = so just a flag
                    name = x
                    val = None
                 else:
                    name = x[:p]
                    val = x[p+1:].strip()
                 if val: 
                    GVars.defines.append("-D%s=%s" % (name.upper(),val))
                 else:
                    GVars.defines.append("-D%s" % name.upper())
                 GVars.definesNames.append("-D%s" % name.upper())
        elif arg == "--particlemethods":
            keywordParticleType = "TYPE"
            keywordOverrideList = ["INIT","MAP","ADV"] 
            nameValuePairs = []
            particleType = ""

            for x in val.split(","):
                x = x.strip()
                p = x.find("=")
                if p < 0: # We require keywords and associated values.
                    raise SetupError("particlemethods option: must contain 'keyword=value' pairs!")
                else:
                    name = x[:p]
                    val = x[p+1:].strip()
                    if (name.upper() == keywordParticleType):
                        particleType = val[:].upper()
                    else:
                        #Only add the name,value pair if it is an option
                        #we can override.
                        if name.upper() in keywordOverrideList:
                            nameValuePairs.append((name.upper(),val.upper()))
                        else:
                            raise SetupError("particlemethods option: keyword '%s' no recognized!" % name)

            if (particleType != ""):
                GVars.particleMethods[particleType] = nameValuePairs
            else:
                raise SetupError("particlemethods option: must specify a particle type with TYPE=<name of particle type>!")

        else:
            if not val: 
               a = arg
            else: a = "%s=%s" % (arg,val)
            GVars.out.put('Invalid Option: %s' % a,globals.ERROR)
            usage()

    GVars.withUnits = withUnits.keys() # takes care of duplicate --with-unit=X arguments

    GVars.out.setDebugLevel(GVars.verbose) # inform GVars.out about the verbosity level

# Perform basic checks on user options 
def checkOpts():
    if GVars.auto and GVars.unitsFile:
        raise SetupError("-unitsfile cannot be used with -auto")

    if ((GVars.gridInterpolation == globals.GRID_INTERP_NATIVE) and
        (GVars.gridGeometry == globals.GRID_GEOM_CYLINDRICAL or
         GVars.gridGeometry == globals.GRID_GEOM_SPHERICAL or
         GVars.gridGeometry == globals.GRID_GEOM_POLAR)):
        raise SetupError("Native Grid interpolation may not be used with %s" % GVars.gridGeometry)
 
    if GVars.dimension < 3: 	 
        if GVars.nzb != None: 	 
            raise SetupError("Must not specify nzb for dimensionality < 3d") 	 
    if GVars.dimension < 2: 	 
        if GVars.nyb != None: 	 
            raise SetupError("Must not specify nyb for dimensionality < 2d") 	 

# finalize all the options made
# this is called after unitList.adjustOpts()
def finalizeOpts():
    if GVars.nzb == None: 	 
      if GVars.dimension > 2: 	 
        GVars.nzb = 8 	 
      else: 	 
        GVars.nzb = 1 	 

    if GVars.nyb == None: 	 
      if GVars.dimension > 1: 	 
        GVars.nyb = 8 	 
      else: 	 
        GVars.nyb = 1 	 
 	 
    if GVars.nxb == None: 	 
       GVars.nxb = 8 	 

    defines = {}
    defines["N_DIM"] = GVars.dimension
    defines["MAXBLOCKS"] = GVars.maxblocks
    defines["NXB"] = GVars.nxb
    defines["NYB"] = GVars.nyb
    defines["NZB"] = GVars.nzb
    
    GVars.setupVars.set("nxb",GVars.nxb)
    GVars.setupVars.set("nyb",GVars.nyb)
    GVars.setupVars.set("nzb",GVars.nzb)
    GVars.setupVars.set("maxBlocks",GVars.maxblocks)

##    if GVars.setupVars.get("ParameshLibraryMode"):
##        defines["LIBRARY"] = None

    ditems = defines.items()
    ditems.sort()
    for (k,v) in ditems:
      if v: 
         GVars.defines.append("-D%s=%s"% (k,v))
      else: GVars.defines.append("-D%s"% k)
      GVars.definesNames.append("-D%s"% k)

    # if old defines file exists, check if new defines conflict with old one
    newDefines = None # the defines have not changed

    od = os.path.join(GVars.flashHomeDir,GVars.objectDir)
    sd = os.path.join(od,globals.SetupDefinesFilename)
    if not os.path.isdir(od): os.makedirs(od)
    f = lazyFile.LazyFile(sd)
    for a in GVars.defines: f.write(a+"\n")
    f.close()
    newDefines = not f.samefile

    # if defines changed and we wanted noClobber
    if newDefines and GVars.noClobber == 1:
       GVars.out.put("Compile Time Parameters changed since last run. Ignoring noclobber",globals.WARN)
       GVars.noClobber = 0

################################## External Interface ################################3

def init():
    pass

def parse():
    parseCommandLine()
    checkOpts()

def final():
    finalizeOpts()

def writeCmdLine():
    file = open(os.path.join(GVars.flashHomeDir,GVars.objectDir,globals.SetupCallFilename), 'w')
    for arg in sys.argv:
        file.write(arg+' ')
    file.write('\n')
    file.write('\nExpanded Command line\n')
    file.write('%s\n' % GVars.fullCmdLine)
    file.write('\nDefined Setup Variables:\n\n')
    GVars.setupVars.printvars(file,ignoreprefix="with")
    file.close()

