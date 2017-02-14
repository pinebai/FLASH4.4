"""Declares global symbols for use in FLASH setup"""

##################################################################
########### Internal Implementation ##############################
##################################################################

import string
import sys
import os.path
import types

# Different levels of verbosity
PPDEBUG = 10  # Used for debugging Config file preprocessing
PPWARN = 20   # Prints preprocessor warnings as well
DEBUG = 30    # Use for "current value of variable x is this
INFO = 40     # Information about what the code is doing (not needed on normal runs)
WARN = 50     # WARNINGs to the user
IMPINFO = 60  # Important information which user should usually be told
ERROR = 70    # Critical Error which must be flagged to user
ALWAYS = 80   # Information which must always be displayed (SUCESS / FAILURE)

# Different grid interpolations
GRID_INTERP_NATIVE = "GRID_INTERP_NATIVE"
GRID_INTERP_MONOTONIC = "GRID_INTERP_MONOTONIC"

# Different geometries
GRID_GEOM_CARTESIAN = "GRID_GEOM_CARTESIAN"
GRID_GEOM_CYLINDRICAL = "GRID_GEOM_CYLINDRICAL"
GRID_GEOM_SPHERICAL = "GRID_GEOM_SPHERICAL"
GRID_GEOM_POLAR = "GRID_GEOM_POLAR"
GRID_GEOM_UNDEF = "GRID_GEOM_UNDEF"

# This string is used many times in the code
# reduces chance of typos
SIM_SETUP_DIR = 'Simulation/SimulationMain'
SETUP_UNITS_FILENAME = 'setup_units'
UNITS_FILENAME = 'Units'
SETUP_PARAMS_FILENAME = "setup_params"
SETUP_VARS_FILENAME = "setup_vars"
SETUP_DEFINES_FILENAME = "setup_defines"
SETUP_CALL_FILENAME = "setup_call"
FLASH_DEFINES_FILENAME = "Flash.h"
SIM_INT_TO_STR_FILENAME = 'Simulation_mapIntToStr.F90'
SIM_STR_TO_INT_FILENAME = 'Simulation_mapStrToInt.F90'
SIM_PARTICLES_VAR_FILENAME = 'Simulation_mapParticlesVar.F90'
BUILD_STAMP_GEN_FILENAME = 'make_bstamp'
SETUP_LIBRARIES_FILENAME = 'setup_libraries'
RP_INIT_PARAMS_FILENAME = 'rp_initParameters.F90'
RP_DEFAULT_PAR_FILENAME = 'default.par'
SETUP_FLASH_UNITS_FILENAME = 'setup_getFlashUnits.F90'
RENORM_GROUP_FILENAME = 'Simulation_getRenormGroup.F90'
VARNAME_TYPE_FILENAME = 'Simulation_getVarnameType.F90'
PM3_RP_FILENAME = 'amr_runtime_parameters'
SUCCESS_FILENAME = ".success"
SIM_PARTICLE_TYPE_FILENAME = 'Particles_specifyMethods.F90'
EOS_MAP_FILENAME = 'eos_variableMap.F90'

# Names of template files. Usually located in bin directory
MAKEFILE_TEMPLATE = 'Makefile.tpl'
BUILD_STAMP_TEMPLATE = 'BuildStamp.tpl'
FLASH_DEFINES_TEMPLATE = 'Flash.h.tpl'
SIM_INT_TO_STR_TEMPLATE = 'MapIntStr.F90.tpl'
SIM_STR_TO_INT_TEMPLATE = 'MapStrInt.F90.tpl'
SIM_PARTICLES_VAR_TEMPLATE = 'getParticlesVarMap.F90.tpl'
SETUP_FLASH_UNITS_TEMPLATE = 'FlashUnits.F90.tpl'
MAKEFILE_STUB_TEMPLATE = 'MakefileStub.tpl'
RENORM_GROUP_TEMPLATE = 'GetRenormGroup.F90.tpl'
VARNAME_TYPE_TEMPLATE = 'GetVarnameType.F90.tpl'
PM3_RP_TEMPLATE = 'amr_runtime_parameters.tpl'
SIM_PARTICLE_TYPE_TEMPLATE = 'Particles_specifyMethods.F90.tpl'
EOS_MAP_TEMPLATE = 'eos_variableMap.F90.tpl'

# List of prefixes for files created by setup. files with following prefixes
# should not be touched by noClobber
NO_CLOBBER_EXCEPTION_LIST = ["setup_",
                             FLASH_DEFINES_FILENAME,
                             SIM_INT_TO_STR_FILENAME,
                             SIM_STR_TO_INT_FILENAME,
                             UNITS_FILENAME,
                             BUILD_STAMP_GEN_FILENAME,
                             RP_INIT_PARAMS_FILENAME,
                             RENORM_GROUP_FILENAME,
                             ]

SHORTCUT_CHAR = "+"

# for processing internal and external libraries
COMPILERS = ('FFLAGS', 'CFLAGS', 'LFLAGS', 'LIB')
DEFLTFLAG = 'OPT'

# unit prefixes which setup knows about
GRID_PREFIX = "Grid/GridMain"
GRID_CHOICES = ["UG", "paramesh", "Samrai", "Chombo"]
SIMULATION_PREFIX = "Simulation/SimulationMain"


class SetupError(Exception):
    """Class for SetupError Exception"""
    pass


class IndentedOutput:
    """Pretty printing class"""
    def __init__(self, num_spaces=4, file=None, debuglevel=WARN):
        self.WRAP = 80
        self.indent = 0
        self.num_spaces = num_spaces
        self.debuglevel = debuglevel

        if not file:
            self.file = sys.stdout
        else:
            self.file = file

    def __pretty(self, text, wrap):
        """Return a list of strings so that each string has length
        <= wrap and consist of full words only. If one word is too
        long it may be an entry in the list."""
        # soln so far
        ans = []

        more = text
        while more:
            if len(more) <= wrap:
                ans.append(more)
                break

            # find last space before or equal to wrap length
            a = more.rfind(" ", 0, wrap + 1)
            if a < 0:
                # first word is long
                b = more.find(" ")
                if b < 0:
                    ans.append(more)
                    more = ""
                else:
                    ans.append(more[:b])
                    more = more[b:]
            else:
                # a >= 0
                ans.append(more[:a])
                more = more[a:]

            # if anything left remove leading spaces
            more = more.strip()

        return ans

    def put(self, text, level=WARN):
        # info not important
        if level < self.debuglevel:
            return

        for line in map(string.strip, string.split(text, '\n')):
            nsp = self.indent * self.num_spaces
            ans = self.__pretty(line, self.WRAP - nsp)

            # printing an empty string
            if not ans:
                self.file.write('\n')
                self.file.flush()
                continue

            for a in ans:
                self.file.write("%s%s\n" % (' ' * nsp, a))
                self.file.flush()

    def pop(self, numLevs=1):
        """Remove numLevs indentation levels"""
        self.indent = self.indent - numLevs
        if self.indent < 0:
            raise SetupError('Popped indentation one time too many!')

    def push(self, numLevs=1):
        self.indent += numLevs


class SetupVarsClass:
    """Handle setup variables"""

    def __init__(self):
        # actually contains the variables and values
        self.vars = {}

        # fbs is true by default
        self.set("fixedBlockSize", True)

    def addunit(self, unitname, val=False):
        """Add a with variable for toplevel unit of specified unit."""
        parts = unitname.split(os.sep)

        while parts and parts[0] and not parts[0][0].isupper():
            del parts[0]

        if parts:
            self.set("with%s" % parts[0], val)

    def set(self, key, value):
        """Set key to value, value may be a quoted string"""
        if not value:
            value = ""

        # direct store
        if type(value) != types.StringType:
            self.vars[key] = value
            return

        # we have a string which may require type change
        if len(value) >= 2 and value[0] == value[-1] and value[0] in ["'", '"']:
            value = str(value[1:-1])
        elif value.lower() in ["true", "false", "yes", "no"]:
            # make it boolean
            value = value.lower() in ["true", "yes"]
        else:
            try:
                a = int(value)
            except:
                pass
            else:
                value = a

        self.vars[key] = value

    def get(self, key):
        return self.vars.get(key, "")

    def printvars(self, file, ignoreprefix="with"):
        """File is a file object opened for writing."""
        vi = self.vars.items()
        vi.sort()
        for k, v in vi:
            if ignoreprefix and k.startswith(ignoreprefix):
                continue
            file.write('%s --> %s %s\n' % (k, v, type(v)))
        file.write('\n')


########################
### Global Variables ###
########################
class GVarsClass(object):
    """Stores Global variables (visible to most of the code). Also includes parsed version of
    command line options"""

    _defaults = {'project_setup_dir': 'setup',
                 'project_build_dir': 'build',
                 'project_source_dir': 'source',
                 'project_simulations_dir': 'simulations',
                 'run_id': None,
                 'run_dir': None,
                 }

    def init(self, flash_src_dir, hard=False):
        self.out = IndentedOutput()        # pretty printer
        self.setup_vars = SetupVarsClass()  # handles setup variables
        self.indexReorder = False          # reorder indices in unk or not

        # set or reset default values
        for k, v in self._defaults.items():
            if not hasattr(self, k) or hard:
                setattr(self, k, v)

        self.project_base_dir = os.path.commonprefix([self.project_source_dir, 
                                                      self.project_simulations_dir])
        self.project_base_dir = os.path.split(self.project_base_dir)[0]
        if self.project_base_dir == '' or not os.path.exists(self.project_base_dir):
            self.project_base_dir = None

        self.portable = 0    # copy or link files
        self.verbose = WARN  # level of verbose printing
        self.report = 0

        self.dimension = 2  # dimensionality of problem trying to solve
        self.defines = []   # PreProcessor symbols to inform compiler
        self.definesNames = [] # PreProcessor symbols to inform compiler, names only
        self.maxblocks = None
        self.makedisplay = 1   # dont display lines as they are executed in makefile
        self.makefileext = ""  # link to Makefile.h+Thisextension
        self.nxb = None  # domain size
        self.nyb = None
        self.nzb = None
        self.build_site = None
        self.build_tau = None
        self.build_os = None
        self.desc_filename = 'flash_desc.json'
        self.auto = 0
        self.simulation_name = None
        self.simulation_dir = None
        self.run_dir_prefix = 'run-'
        self.buildFlag = 'OPT'  # other options are 'DEBUG' and 'TEST' and 'OPT'
        self.npg = 0
        self.withUnits = []
        self.with_libraries = {}
        self.noClobber = 0

        # list of unix patterns for files to be copied over to object directory
        self.datafiles = []

        # list of unix patterns for files to be copied over to object directory
        self.parfile = ""

        # units to be removed after all conditions have been met (USER BEWARE)
        self.killUnits = {}
        self.withoutUnits = {}  # units not needed (unless really needed)
        self.unitsfile = ""  # name of Units file to pick the units from

        # for completeness store all withoutLibrary options (now dont do anyting with it)
        self.withoutLibraries = {}

        self.gridInterpolation = None
        self.gridGeometry = GRID_GEOM_UNDEF
        self.curvilinear = 0
        self.strictParams = 0
        self.particleMethods = {}
        self.particlesUnitIncluded = False
        self.nonexistent = "NONEXISTENT"

        self.flash_src_dir = flash_src_dir
        self.source_dir = os.path.join(self.flash_src_dir, 'source')
        self.source_path = [self.project_source_dir, self.source_dir]
        self.simulations_dir = os.path.join(self.source_dir, SIM_SETUP_DIR)
        self.simulations_path = [self.project_simulations_dir] + \
                                [os.path.join(sp, SIM_SETUP_DIR) for sp in self.source_path]

        self.libDir = os.path.join(self.flash_src_dir, 'lib')
        self.binDir = os.path.join(self.flash_src_dir, 'bin')
        self.tauDir = os.path.join(self.flash_src_dir, 'tools', 'tau')

        self.topUnitNames = self.get_top_unit_names()
        self.unit_aliases = {self.project_simulations_dir: SIM_SETUP_DIR}

        self.extra_runtime_files = ['amr_runtime_parameters']

    def get_top_unit_names(self, path=None, ignoreUnits=("flashUtilities",), init_path=None):
        """Return the list of all top level unit names.
        path may be a list of str, a str, or None.  When None (default)
        the value from self.source_path is taken.
        """
        if path is None:
            path = self.source_path
        elif isinstance(path, basestring):
            path = [path]

        if init_path is None:
            init_path = path

        # verify that we have at least one directory to look in
        if not any([os.path.isdir(p) for p in path]):
            return []

        # Find units 
        ans = []
        for p in path:
            if not os.path.exists(p):
                continue

            basename = os.path.basename(os.path.abspath(p))
            ch = basename[0]

            # Are we a unit
            if 'A' <= ch <= 'Z':
                return [os.path.relpath(p, ip) for ip in init_path if p.startswith(ip)]

            # we are only an organizational directory
            for dir in os.listdir(p):
                if dir[0] == "." or (dir in ignoreUnits):
                    continue

                subdir = os.path.join(p, dir)
                if not os.path.isdir(subdir):
                    continue

                ans.extend(self.get_top_unit_names(subdir, ignoreUnits, init_path))

        # make unique and return
        ans = sorted(set(ans))
        return ans


#############################################################################
##################### initialization code   #################################
#############################################################################

gvars = GVarsClass()
