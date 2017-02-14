"""Code which parses the command line arguments for the setup command."""

import string
import re
import getopt
import sys
import os
import shlex

# Local imports
# Relative imports needed
from .. import FLASH_SRC_DIR
from . import setup_globals
from .lazy_file import LazyFile
from .setup_globals import gvars, SetupError
from .utils import dir_glob, search_paths

if sys.version_info[0] == 3 or sys.version_info[:2] == (2, 7):
    import argparse
else:
    from .. import _argparse as argparse

from .utils import add_simple_opts

WHITESPACE = re.compile('\s')


########################## Internal implementation ##########################

# default arguments
ADDL_DEF_ARGS = ["+default"]

# options without arguments
WITHOUT_ARGS = ["auto","1d","2d","3d","portable",
                "makehide", "curvilinear",
                "opt","debug","test", "index-reorder", "strictparams",
                "fbs","help", "noclobber", "nofbs"] # nofbs disables fixedBlockSize
# options with arguments
WITH_ARGS = ["maxblocks","nxb","nyb","nzb","verbose","site","ostype",
             "defines", "objdir", "locsrcdir", "rundir", "with-unit", "unit", "with-library",
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
       if not hasattr(gvars, 'shortcuts'):
           gvars.init(FLASH_SRC_DIR)
           gvars.shortcuts = getShortcuts()
       sitems = gvars.shortcuts.items()
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
   if hasattr(gvars, 'out'):
        gvars.out.put("Processing Shortcut file: %s" % filename,setup_globals.IMPINFO)
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
      if hasattr(gvars, 'out'):
         gvars.out.put("No Shortcut file specified using default",setup_globals.DEBUG)
      sfiles = [FLASH_SRC_DIR + "/bin/setup_shortcuts.txt"]
   else:
      sfiles = os.environ["SETUP_SHORTCUTS"].split(":")
   shortcuts = {}
   for name in sfiles:
      sfile = os.path.abspath(name)
      if not os.path.exists(sfile):
         if hasattr(gvars, 'out'):
            gvars.out.put("Unable to open %s. Ignoring" % sfile,setup_globals.WARN)
      else: shortcuts.update(getShortcutDict(sfile))
   return shortcuts

def expandShortcuts(args,RecLimit=255):
   """Expand all shortcuts specified in args. 
      No more than RecLimit expansions should be required. 
      This limit can be used to detect circular references"""
   ans = []
   count = RecLimit
   gvars.shortcuts = getShortcuts()
   shortcuts = gvars.shortcuts
   while args:
      scut = args[0]
      del args[0]
      if scut[0] != setup_globals.SHORTCUT_CHAR: 
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
            gvars.out.put("\n***WARNING*** Ignoring unknown shortcut %s while expanding %s.\n" % (scut," ".join(args)),setup_globals.IMPINFO)
         elif count > 0:
            count = count - 1
            args[0:0] = shortcuts[scut] # insert expansion of shortcut in front of args
         else:
            raise SetupError("Too many shortcuts. Circular reference?")
   return ans

def setSetupVars(strlist):
   """Given a list of strings. Identifies those of the form A=B and 
     updates gvars.setup_vars.  Returns unprocessed list of strings"""
   ans = []
   # declare setup variables for all top level units
   # these are boolean variables which decide whether
   # to include this unit or not
   for x in gvars.topUnitNames: # by default do not include any units
       gvars.setup_vars.addunit(x,False)
   for opt in strlist:
       p = opt.find("=")
       if p < 0: # no equal found
          ans.append(opt)
       else:
          gvars.setup_vars.set(opt[:p],opt[p+1:]) # set variable to value
   return ans
 
def custom_getopt(args, longoptions):
   """Parse command line arguments, working in gnu mode using only longoptions.
     GNU MODE =  allow non-option arguments to be followed by options.
     We also allow user to use shortcuts as well as define variables like make"""
   max_shortcuts = 255 # no more than 255 shortcuts allowed per invocation
   optvallist = []
   rest = []
   toparse = cleanupCommandLine(ADDL_DEF_ARGS+args) # double "-" to "--"
   fulllist = expandShortcuts(toparse) # expand the shortcuts
   fullCmdList = [x for x in fulllist if x and (x[0] == "-" or x.find("=") < 1) ] # remove non-options with "=" in it
   gvars.fullCmdLine = " ".join(fullCmdList) # this gets printed
   toparse = fulllist[:] # make a copy of the list
   while toparse:
    optval_sublist, toparse = getopt.getopt(toparse,"",longoptions)
    optvallist.extend(optval_sublist)
    if toparse: 
       badoption = toparse[0]
       del toparse[0]
       rest.append(badoption)
   return (fulllist,optvallist,setSetupVars(rest))

# FIXME So many bugs in trying to merge setup arguments. 
# Forgetting about it for the moment
#def merge_args(arg_lists):
#    """Merge the argument lists for setup."""
#    arg_lists = [(shlex.split(a) if isinstance(a, basestring) else a) for a in arg_lists]
#    # generate parser for setup
#    setup_parser = argparse.ArgumentParser(prefix_chars="+-")
#    for arg_list in arg_lists:
#        add_simple_opts(setup_parser, arg_list)
#    setup_parser.add_argument('positionals', type=str, nargs=argparse.REMAINDER)
#
#    # parse the options into the namespace, order matters!
#    ns = argparse.Namespace()
#    for arg_list in arg_lists:
#        ns = setup_parser.parse_args(arg_list, ns)
#
#    # create the qsub command line string from the namespace
#    opt_strings = {act.dest: act.option_strings[0] for act in setup_parser._actions}
#    merged_args = ['qsub']
#    for opt, value in ns.__dict__.items():
#        merged_args.append(opt_strings[opt])
#        if value is None:
#            continue
#        elif WHITESPACE.search(value):
#            value = '"' + value + '"'
#        merged_args.append(value)
#    return merged_args

def merge_args(arg_lists):
    margs = [(shlex.split(a) if isinstance(a, basestring) else a) for a in arg_lists if 0 < len(a)]
    if 0 == len(margs):
        raise ValueError("no setup arguments given.")
    return margs[0]

def parseCommandLine(arg_lists):
    # process all the options
    # those requiring arguments suffixed with =
    longopts = WITHOUT_ARGS + [x+"=" for x in WITH_ARGS]

    merged_args = merge_args(arg_lists)

    try:
       (fullcmdline,optvallist,rest) = custom_getopt(merged_args, longopts)
    except getopt.GetoptError,e:
       gvars.out.put(str(e),setup_globals.ERROR)
       usage() # print usage info

    # Not given any program name
    if not rest:
        usage()

    # if given many the last one wins
    if len(rest) > 1:
       gvars.out.put("         WARNING: Multiple problem names given. winner="+rest[-1],setup_globals.WARN)
    gvars.simulation_name = rest[-1] # store program name

    # Find the approrpiate simualtion directory 
    # a little tricky to enable project_simulations_dir
    sim_dir = search_paths(gvars.simulations_path, gvars.simulation_name)
    if not sim_dir.startswith(gvars.project_simulations_dir):
        sim_dir = [os.path.relpath(sim_dir, sp) for sp in gvars.source_path 
                                                if sim_dir.startswith(sp)][0]
    gvars.simulation_dir = sim_dir 

    # all acceptable values for "--verbose" keyword
    vrblevels = {"DEBUG":setup_globals.DEBUG,
                 "PPDEBUG":setup_globals.PPDEBUG,
                 "PPWARN":setup_globals.PPWARN,
                 "WARN":setup_globals.WARN,
                 "INFO":setup_globals.INFO,
                 "IMPINFO":setup_globals.IMPINFO,
                 "ERROR":setup_globals.ERROR}

    # all acceptable values for "--gridinterpolation" keyword
    allGridInterpolations = {"MONOTONIC":setup_globals.GRID_INTERP_MONOTONIC,
                             "NATIVE":setup_globals.GRID_INTERP_NATIVE}

    # all acceptable values for "--geometry" keyword
    allGeometries = {"CARTESIAN":setup_globals.GRID_GEOM_CARTESIAN,
                     "CYLINDRICAL":setup_globals.GRID_GEOM_CYLINDRICAL,
                     "SPHERICAL":setup_globals.GRID_GEOM_SPHERICAL,
                     "POLAR":setup_globals.GRID_GEOM_POLAR}

    withUnits = {} # dictionary to store list of units to add (prevents duplicate entries)

    for (arg,val) in optvallist:
        if   arg == '--portable':         gvars.portable = 1
        elif arg == '--noclobber':        gvars.noClobber = 1
        elif arg == '--auto':             gvars.auto = 1
        elif arg == '--datafiles':        gvars.datafiles.append(val)
        elif arg == '--parfile':          gvars.parfile = val
        elif arg == '--debug':            gvars.buildFlag = "DEBUG"
        elif arg == '--test':             gvars.buildFlag = "TEST"
        elif arg == '--opt':              gvars.buildFlag = "OPT"
        elif arg == '--maxblocks':          gvars.maxblocks = int(val)
        elif arg == '--makehide':           gvars.makedisplay = 0
        elif arg == '--makefile':           gvars.makefileext = "."+val
        elif arg == '--nxb':                gvars.nxb = int(val)
        elif arg == '--nyb':                gvars.nyb = int(val)
        elif arg == '--nzb':                gvars.nzb = int(val)
        elif arg == '--site':               gvars.build_site = val
        elif arg == '--tau':                gvars.build_tau = val
        elif arg == '--ostype':             gvars.build_os = val
        elif arg == '--objdir':             gvars.project_build_dir = val
        elif arg == '--rundir':             gvars.rundir = val
        elif arg == '--help':               usage()
        elif arg == '--unitsfile':          gvars.unitsfile = val
        elif arg == '--fbs':                gvars.setup_vars.set("fixedBlockSize",True)
        elif arg == '--nofbs':              gvars.setup_vars.set("fixedBlockSize", False)
        elif arg == '--strictparams':       gvars.strictParams = 1
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
        # deprecated with introduction of 'gridinterpolation' (below)
        # DEV 'curvilinear' now deprecated with introduction of 'gridinterpolation' (below)
        elif arg == '--curvilinear':
             gvars.out.put("\n***************************** INFO *********************************\n" +
                           "The -curvilinear flag is nearly always unnecessary.\n" +
                           "Use --gridinterpolation={cylindrical,spherical,polar} to configure\n" +
                           "FLASH for a specific non-Cartesian geometry.\n" +
                           "(These flags are not needed at all when using PARAMESH in LIBRARY\n" +
                           "MODE, i.e., ParameshLibraryMode=True, or when using Paramesh4dev.)\n" +
                           "********************************************************************\n",setup_globals.WARN)
             gvars.curvilinear = 1
        elif arg in ["--1d","--2d","--3d"]: 
             gvars.dimension = int(arg[2]) # i.e. 1 or 2 or 3
             gvars.setup_vars.set("nDim",gvars.dimension)
        elif arg == "--index-reorder":            
             gvars.setup_vars.set("GridIndexReordered",True) # for use in Config file
             gvars.indexReorder = True # for use in Makefile
             gvars.defines.append("-DINDEXREORDER") # for use in code
             gvars.definesNames.append("-DINDEXREORDER") # for use in code
        elif arg in ["--with-unit","--unit"]:
             if val.endswith(os.sep): val = val[:-1]
             if val.startswith("source"+os.sep): val = val[7:]
             withUnits[val] = 1
             gvars.setup_vars.addunit(val,True) # set corresponding setup var to true
        elif arg == "--without-unit": # remove all units added which come under specified unit
             # also add to list of units to be ignored when handling REQUESTS keyword
             if val.endswith(os.sep): val = val[:-1] 
             rmlist = [x for x in withUnits.keys() if x == val or x.startswith(val+os.sep)]
             for x in rmlist:
                 del withUnits[x]
             gvars.withoutUnits[val] = 1
        elif arg == "--kill-unit":
             gvars.killUnits[val] = 1
        elif arg in ["--with-library","--library"]: 
             if not val: continue # no name given --> ignore 
             parts = val.split(",")
             libname = parts[0]
             args = string.join(parts[1:]," ") # replace commas with space
             if len(args) >= 2 and args[0] == args[-1] and args[0] in ['"',"'"]: # argument has been quoted
                args = args[1:-1]
             gvars.with_libraries[libname.lower()] = args
        elif arg == "--without-library":
             val = val.lower()
             if gvars.with_libraries.has_key(val): del gvars.with_libraries[val]
             gvars.withoutLibraries[val] = 1
        elif arg == '--verbose': # set verbosity level
             if not val: continue # no argument dont change level
             if vrblevels.has_key(val.upper()): 
                gvars.verbose = vrblevels[val.upper()]
             else:
                gvars.out.put("Unrecognized verbosity level [%s]" % val,setup_globals.ERROR)
                usage()
        elif arg == '--gridinterpolation':  # set grid interpolation
             if not val: continue  # no argument; don't change interpolation
             if allGridInterpolations.has_key(val.upper()):
                 gvars.gridInterpolation = allGridInterpolations[val.upper()]
             else:
                 gvars.out.put("Unrecognized grid interpolation [%s]" % val, setup_globals.ERROR)
                 usage()                 
        elif arg == '--geometry': # set geometry
             if not val: continue # no argument; don't change geom
             if allGeometries.has_key(val.upper()):
                 gvars.gridGeometry = allGeometries[val.upper()]
                 if val.upper() != "CARTESIAN":
                     # All geometries other than cartesian (i.e. cylindrical,
                     # spherical, and polar) are automatically curvilinear.
                     # (See note above on -curvilinear.)
                     gvars.curvilinear = 1
             else:
                 gvars.out.put("Unrecognized geometry [%s]" % val, setup_globals.ERROR)
                 usage()
        elif arg == "--defines": # declare additional CPP/FPP stuff
             if not val: # kill existing defines
                gvars.defines = []
                gvars.definesNames = []
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
                    gvars.defines.append("-D%s=%s" % (name.upper(),val))
                 else:
                    gvars.defines.append("-D%s" % name.upper())
                 gvars.definesNames.append("-D%s" % name.upper())
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
                gvars.particleMethods[particleType] = nameValuePairs
            else:
                raise SetupError("particlemethods option: must specify a particle type with TYPE=<name of particle type>!")

        else:
            if not val: 
               a = arg
            else: a = "%s=%s" % (arg,val)
            gvars.out.put('Invalid Option: %s' % a,setup_globals.ERROR)
            usage()

    # takes care of duplicate --with-unit=X arguments
    gvars.withUnits = withUnits.keys() 

    # inform gvars.out about the verbosity level
    gvars.out.debuglevel = gvars.verbose

# Perform basic checks on user options 
def checkOpts():
    if gvars.auto and gvars.unitsfile:
        raise SetupError("-unitsfile cannot be used with -auto")

    if ((gvars.gridInterpolation == setup_globals.GRID_INTERP_NATIVE) and
        (gvars.gridGeometry == setup_globals.GRID_GEOM_CYLINDRICAL or
         gvars.gridGeometry == setup_globals.GRID_GEOM_SPHERICAL or
         gvars.gridGeometry == setup_globals.GRID_GEOM_POLAR)):
        raise SetupError("Native Grid interpolation may not be used with %s" % gvars.gridGeometry)
 
    if gvars.dimension < 3: 	 
        if gvars.nzb != None: 	 
            raise SetupError("Must not specify nzb for dimensionality < 3d") 	 
    if gvars.dimension < 2: 	 
        if gvars.nyb != None: 	 
            raise SetupError("Must not specify nyb for dimensionality < 2d") 	 

# finalize all the options made
# this is called after unit_list.adjust_opts()
def finalizeOpts():
    if gvars.nzb == None: 	 
      if gvars.dimension > 2: 	 
        gvars.nzb = 8 	 
      else: 	 
        gvars.nzb = 1 	 

    if gvars.nyb == None: 	 
      if gvars.dimension > 1: 	 
        gvars.nyb = 8 	 
      else: 	 
        gvars.nyb = 1 	 
 	 
    if gvars.nxb == None: 	 
       gvars.nxb = 8 	 

    defines = {}
    defines["N_DIM"] = gvars.dimension
    defines["MAXBLOCKS"] = gvars.maxblocks
    defines["NXB"] = gvars.nxb
    defines["NYB"] = gvars.nyb
    defines["NZB"] = gvars.nzb
    
    gvars.setup_vars.set("nxb",gvars.nxb)
    gvars.setup_vars.set("nyb",gvars.nyb)
    gvars.setup_vars.set("nzb",gvars.nzb)
    gvars.setup_vars.set("maxBlocks",gvars.maxblocks)

##    if gvars.setup_vars.get("ParameshLibraryMode"):
##        defines["LIBRARY"] = None

    ditems = defines.items()
    ditems.sort()
    for (k,v) in ditems:
      if v: 
         gvars.defines.append("-D%s=%s"% (k,v))
      else: gvars.defines.append("-D%s"% k)
      gvars.definesNames.append("-D%s"% k)

    # if old defines file exists, check if new defines conflict with old one
    newDefines = None # the defines have not changed

    #od = os.path.join(gvars.flash_src_dir,gvars.project_build_dir)
    od = os.path.join(gvars.flash_src_dir, gvars.project_setup_dir)
    sd = os.path.join(od,setup_globals.SETUP_DEFINES_FILENAME)
    if not os.path.isdir(od): os.makedirs(od)
    f = LazyFile(sd)
    for a in gvars.defines: f.write(a+"\n")
    f.close()
    newDefines = not f.samefile

    # if defines changed and we wanted noClobber
    if newDefines and gvars.noClobber == 1:
       gvars.out.put("Compile Time Parameters changed since last run. Ignoring noclobber",setup_globals.WARN)
       gvars.noClobber = 0

################################## External Interface ################################3

def parse(arg_lists):
    parseCommandLine(arg_lists)
    checkOpts()

def final():
    finalizeOpts()

def write_cmd_line():
    #file = open(os.path.join(gvars.flash_src_dir, gvars.project_build_dir, setup_globals.SETUP_CALL_FILENAME), 'w')
    file = open(os.path.join(gvars.project_setup_dir, setup_globals.SETUP_CALL_FILENAME), 'w')
    for arg in sys.argv:
        file.write(arg+' ')
    file.write('\n')
    file.write('\nExpanded Command line\n')
    file.write('%s\n' % gvars.fullCmdLine)
    file.write('\nDefined Setup Variables:\n\n')
    gvars.setup_vars.printvars(file,ignoreprefix="with")
    file.close()

