"""This module contains code to pre-process text files, code for a FalshUnit class
which parses Config files, code for FlashUnitUnion class to merge FlashUnits, and
a FlashUnitList class which contains a list of unit names together with methods to
manipulate them."""

import re
import os
import sys
import copy
import collections
from warnings import warn

if sys.version_info[0] == 3 or sys.version_info[:2] == (2, 7):
    from argparse import Namespace
else:
    from .._argparse import Namespace

# relative imports needed
from ..utils import strip_comments

#
# Configuration code
#
RIF = re.compile(r"^\s*IF\s*(?P<cond>.*)\s*$")
RELIF = re.compile(r"^\s*(?:ELSEIF|ELIF)\s*(?P<cond>.*)\s*$")
RANYIF = re.compile(r"^\s*(?:IF|ELSEIF|ELIF)\s*(?P<cond>.*)\s*$")
RERR = re.compile(r"^\s*SETUPERROR\s*(?P<msg>.*)\s*$")
RELSE = re.compile(r"^\s*ELSE\s*(?:#.*)?$")
RENDIF = re.compile(r"^\s*END\s*IF\s*(?:#.*)?$")
RUSE = re.compile(r"^\s*USESETUPVARS\s+(?P<vars>.*)(?:#.*)?$")
RD = re.compile(r"^\s*D\s+")
RPARAM = re.compile(r"^\s*PARAMETER\s+")

# Must register the restriction on keywords
TOP_UNIT_PARSER_SET = set(['D', 'DEFAULT', 'PARAMETER', 'CHILDORDER'])
TOP_UNIT_ERROR = ('Illegal keyword "{0}" in file "{1}", line {2}.\n'
                  "Only DEFAULT, PARAMETER, and 'D' keywords are "
                  "allowed in a top-level Config.")

NO_RESTRICT_PARSER_SET = set([ 'CHILDORDER', 'DATAFILES', 'D', 'FACEVAR', 
    'VARIABLE', 'LIBRARY', 'GUARDCELLS', 'FLUX', 'SCRATCHVAR', 
    'SCRATCHCENTERVAR', 'SCRATCHFACEXVAR', 'SCRATCHFACEYVAR', 'SCRATCHFACEZVAR',
    'SPECIES', 'MASS_SCALAR', 'PARTICLEPROP', 'PARTICLETYPE', 'PARTICLEMAP', 
    'SUGGEST', 'REQUIRES', 'REQUESTS', 'PARAMETER', 'DEFAULT', 'EXCLUSIVE', 
    'KERNEL', 'LINKIF', 'CONFLICTS', 'PPDEFINE', 'NONREP'])
NO_RESTRICT_ERROR = 'Unknown keyword "{0}" in file "{1}", line {2}.'

EOS_STATIC_LIST = ['PRES',  'DENS',  'EINT',  'TEMP',  'GAMC',  'GAME',  'ENER', 
                   'VELX',  'VELY',  'VELZ',  'SUMY',  'YE',    'ENTR',  'PRES1',
                   'PRES2', 'PRES3', 'EINT1', 'EINT2', 'EINT3', 'TEMP1', 'TEMP2',
                   'TEMP3', 'E1',    'E2',    'E3',    'SELE',  'SRAD',]
EOS_STR = "|".join([x.lower() + "|" + x.upper() for x in EOS_STATIC_LIST])


# Config regular expressions
RE_EOS_STR = '(\s+EOSMAP:\s*(?P<eosmap>' + EOS_STR + ')\s*$)?'
RE_EOS_MAP = 'EOSMAP(INOUT)?:\s*(?P<eosmap>' + EOS_STR + \
             ')|EOSMAPIN:\s*(?P<eosmapin>' + EOS_STR + \
             ')|EOSMAPOUT:\s*(?P<eosmapout>' + EOS_STR + ')'
RE_TYPE = 'TYPE:\s*(?P<type>GENERIC|PER_VOLUME|PER_MASS)'

RE_CHILDORDER = re.compile('CHILDORDER\s+(.*)$')
RE_DATAFILES = re.compile('DATAFILES\s+(.*)$')
RE_D = re.compile('D\s*(&|\S+)\s*(.*)$')
RE_FACEVAR = re.compile('FACEVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_VARIABLE = re.compile('VARIABLE\s+(?P<varname>\w+)(\s+)?(?:(' + RE_TYPE + \
                         '|' + RE_EOS_MAP + r')\s*)*\s*$')
RE_LIBRARY = re.compile('LIBRARY\s+(\S+)\s*(.*)$')
RE_GUARDCELLS = re.compile(r'GUARDCELLS\s+([0-9]+)$')
RE_FLUX = re.compile(r'FLUX\s+(\w+)$')
RE_SCRATCHVAR = re.compile('SCRATCHVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_SCRATCHCENTERVAR = re.compile('SCRATCHCENTERVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_SCRATCHFACEXVAR = re.compile('SCRATCHFACEXVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_SCRATCHFACEYVAR = re.compile('SCRATCHFACEYVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_SCRATCHFACEZVAR = re.compile('SCRATCHFACEZVAR\s+(?P<varname>\w+)' + RE_EOS_STR)
RE_SPECIES = re.compile(r'SPECIES\s+(?P<name>\w+)(\s+TO\s+(?P<nElec>\d+))?\s*$')
RE_MASS_SCALAR = re.compile(r'MASS_SCALAR\s+(?P<msname>\w+)(\s+)?(?:(NORENORM|RENORM:'
                            r'(?P<gpname>\S*)|' + RE_EOS_MAP + ')\s*)*\s*$')
RE_PARTICLEPROP = re.compile(r'PARTICLEPROP\s+(?P<name>\w+)\s+(?P<type>INTEGER|REAL)\s*$')
RE_PARTICLETYPE = re.compile(r'PARTICLETYPE\s+(?P<particleType>\w+)\s+INITMETHOD\s+'
                             r'(?P<initMethod>\w+)\s+MAPMETHOD\s+(?P<mapMethod>\w+)'
                             r'\s+ADVMETHOD\s+(?P<advMethod>\w+)\s*$')
RE_PARTICLEMAP = re.compile(r'PARTICLEMAP\s+(TO)\s+(?P<name>\w+)\s+(FROM)\s+(?P<type>'
                            r'(SCRATCHVAR|SCRATCHCENTERVAR|SCRATCHFACEXVAR|SCRATCHFACEYVAR|'
                            r'SCRATCHFACEZVAR|FACEX|FACEY|FACEZ|VARIABLE|MASS_SCALAR|SPECIES'
                            r'))\s+(?P<varname>\w+)\s*$')
RE_SUGGEST = re.compile(r'SUGGEST\s+(.*)$')
RE_REQUIRES = re.compile(r'REQUIRES\s+(?P<orList>(?:\S+|\s+OR\s+)+)\s*$')
RE_REQUESTS = re.compile(r'REQUESTS\s+(?P<name>\S+)\s*$')
RE_PARAMETER = re.compile(r'PARAMETER\s+(\w+)\s+(INTEGER|REAL|STRING|BOOLEAN)'
                          r'\s+(?:(CONSTANT)\s+)?([^[]+)(?:[[](.*)[]])?\s*$')
RE_DEFAULT = re.compile(r'DEFAULT\s+(\S+)$')
RE_EXCLUSIVE = re.compile( r'EXCLUSIVE\s+(.+)$')
RE_KERNEL = re.compile(r'KERNEL(\s+([^ \t]+))?\s*$')
RE_LINKIF = re.compile(r'LINKIF\s+(.+)$')
RE_CONFLICTS = re.compile(r'CONFLICTS\s+(.+)$')
RE_PPDEFINE = re.compile('PPDEFINE\s+(?P<sym>\w+)(?:\s+(?P<val>\w+))?$')
RE_NONREP = re.compile('NONREP\s+(?P<type>MASS_SCALAR|VARIABLE)\s+(?P<name>\w+)'
                       '\s+(?P<nloc>[0-9]+)\s+(?P<pglob>\w+)(\s+(?P<namef>\S+))?$')
# HACK for STRING PARAMETERS since "parse" doesn't handle string parameters
RE_PSTRING = re.compile(r'PARAMETER\s+(\w+)\s+STRING\s+(".*")\s*(?:[[](.*)[]])?(?:#.*)?$')



class Configuration(object):
    """Encapsulates FLASH Unit information as expressed in Config files."""

    def __init__(self, conffile=None, top_unit=False, usevars=None, ignore_pp=False, 
                 unitname="", unitpath='', without_units=None):
        """Data is accessed through the appropriate attributes::

            Configuration('IO/common/Config').default
            Configuration('../Simulation/SimulationMain/Sedov/Config').library
            Configuration('').flux

        Parameters
        ----------
        conffile : str or file object or None
            path to the config file, None for clean instance.
        top_unit : bool, optional
            whether the configuration file is for a top-level FLASH unit.
        usevars : dict or None, optional
            dictionary of initial values for attributes.  Often given 
            as flmake.gvars.setup_vars.vars for setting up actual problems.
        ignore_pp : bool, optional
            Igonore pre-processing for documentation purposes.
        unitname : str, optional
            Name of the unit, if applicable
        unitpath : str, optional
            Absolute path to the unit, if applicable.  Useful for when conffile
            doesn't exist.
        without_units : dict or None, optional
            Units to not allow, in the case of requests.  Bool valued.

        """
        self._init_usevars = {} if usevars is None else usevars
        self.usevars = {}
        self.conffile = conffile
        self.filename = conffile if isinstance(conffile, basestring) else ""
        self.top_unit = top_unit
        self.unitname = unitname
        self.unitpath = unitpath
        self.without_units = {} if without_units is None else without_units

        self.lineno = 0
        self.ignore_pp = ignore_pp
        self.normalkeywords = ["SETUPERROR"]  # keywords to be ignored inside a block

        # initial data vaules
        self.childorder = []
        self.datafiles = []
        self.d = {}
        self.facevar = {}
        self.variable = {}
        self.library = {}
        self.guardcells = None
        self.flux = {}
        self.scratchvar = {}
        self.scratchcentervar = {}
        self.scratchfacexvar = {}
        self.scratchfaceyvar = {}
        self.scratchfacezvar = {}
        self.species = {}
        self.mass_scalar = {}
        self.particleprop = {}
        self.particletype = {}
        self.particlemap = {}
        self.suggest = []
        self.requires = []
        self.requests = []
        self.parameter = {}
        self.default = ""
        self.exclusive = []
        self.kernel = None
        self.linkif = []
        self.conflicts = {}
        self.ppdefine = {}
        self.nonrep = {}

        for key, value in self._init_usevars.items():
            key = key.lower()
            attr = getattr(self, key, None)
            if attr is None:
                setattr(self, key, copy.copy(value))

        # Parse, if we can
        if self.conffile is not None:
            self.parse()

    @property
    def parser_set(self):
        return TOP_UNIT_PARSER_SET if self.top_unit else NO_RESTRICT_PARSER_SET

    @property            
    def error_message(self):
        return TOP_UNIT_ERROR if self.top_unit else NO_RESTRICT_ERROR

    def process_file(self):
        """Processes the config file."""
        opened_here = False
        if isinstance(self.conffile, basestring):
            cf = open(self.conffile, "r")
            opened_here = True
        else:
            cf = self.conffile

        # If the first line of the file starts with ##python',
        # then we compile it as python and run its genLines method to get the lines
        if cf.readline().startswith("##python:genLines"):
            # read the rest of the file (doesnt include ##python:genLines line)
            # and compile it into code
            code = compile(cf.read(), self.filename, "exec")
            cf.close()

            # initialize a namespace with the prelude code
            ns = Namespace() 

            # execute the code object within the prelude'd namespace
            exec code in ns.__dict__

            # helper to split multi-line strings and flatten nested
            # iterables into one long iterator of strings
            def flatten(it):
                if type(it) == type(""):
                    for x in it.split('\n'):
                        yield x
                else:
                    for x in it:
                        for y in flatten(x):
                            yield y

            lines = flatten(ns.genLines(self._init_usevars.copy()))
        else:
            # just a regular config file
            cf.seek(0)
            lines = cf.readlines()

        # exception goes to caller
        if self.ignore_pp:
            # In this case, this method has been called expressly for
            # the purpose of generating documentation, not because we
            # are really setting up a problem, so we use "filter_lines"
            # to strip everything from the Config file but parameter
            # names and values (prefixed by the keyword "PARAMETER")
            # and parameter descriptions (prefixed by a capital 'D')
            ans = self._filter_lines(lines)
        else:
            ans = self._process_lines(lines)

        if opened_here:
            cf.close()
        return ans

    def _pp_keyword(self, line):
        """Returns "IF,ELSE,ENDIF,USE,None" depending on syntax of line."""
        if RIF.match(line):
            return "IF"
        elif RELIF.match(line):
            return "ELIF"
        elif RELSE.match(line):
            return "ELSE"
        elif RENDIF.match(line):
            return "ENDIF"
        elif RUSE.match(line):
            return "USESETUPVARS"
        elif RERR.match(line):
            return "SETUPERROR"
        else:
            return None

    def _check_use(self, line):
        """Confirm that variables named in line are present. Otherwise raise a warning."""
        m = RUSE.match(line)
        if not m:
            return
        vars = [x.strip() for x in m.group("vars").split(",")]
        badlist = [x for x in vars if x not in self.usevars]

        # all declared variables (except for the with* ones)
        initvars = [x for x in self.usevars.keys() if not x.startswith("with")]

        # found missing variables
        if 0 < len(badlist):
            msg = ["\nConfiguration Warning: \nRequired uninitialized variables"
                   " will be initialized to ''"]
            msg.append("File: {0}\nUnnitialized Variables: {1}".format(self.filename,
                                                                       ", ".join(badlist)))
            msg.append("Initialized Variables: {0}".format(", ".join(initvars)))
            msg.append("NOTE: * Variable Names are case sensitive. ")
            msg.append("      * The withUNIT variables are not shown in the initialized list.")
            msg = "\n".join(msg)
            warn(msg, SyntaxWarning)
            for x in badlist:
                self.usevars[x] = ""

    def _handle_pp_cmd(self, cmd, line):
        if cmd == "SETUPERROR":
            # need to raise an error
            merr = RERR.match(line)
            if not merr:
                return
            halt = "Setup halted by rule defined in {0} (line {1}): "
            halt = halt.format(self.filename, self.lineno)
            raise SyntaxError(halt + merr.group("msg"))
        else:
            raise SyntaxError("Programming error. Got unhandled configuration command: " + cmd)

    def _eval_cond(self, line):
        """Evaluates condition and returns True/False/None.
        Returns None if line has bad syntax"""
        
        # either IFDEF or ELSEIF
        m = RANYIF.match(line)
        if not m:
            return None

        cond = m.group("cond")
        valcopy = dict(self.usevars)

        try:
            value = eval(cond, {}, valcopy)
        except Exception as e:
            msg = 'File: {0}\nLine: {1}\nCondition: {2}\nVariables: {3}\nDetails: {4}'
            msg = msg.format(self.filename, self.lineno, cond, self.usevars, e)
            raise SyntaxError(msg)
        else:
            # Warn user about side effect
            if valcopy != self.usevars:
                msg = ["Configuration SIDEEFFECT Warning: Condition has side effect"]
                msg_line = "File: {0}\nLine : {1}\nCode: {2}\nBefore: {3}\nAfter: {4}"
                msg_line = msg_line.format(self.filename, self.lineno, cond, self.usevars, valcopy)
                msg.append(msg_line)
                msg.append("Trashing Side effects")
                msg = "\n".join(msg)
                warn(msg, SyntaxWarning)

        return value

    def _filter_lines(self, lines):
        """Remove everything except for lines starting with 'D '
        (param descriptions) and the keyword 'PARAMETER'"""
        # lineno being on self prevents this from beinf
        ans = []
        for self.lineno, line in enumerate(lines, 1):
            if RD.match(line) or RPARAM.match(line):
                ans.append((self.lineno, line))
        return ans

    def _process_lines(self, lines):
        """Here we evaluate arbitrarily nested conditionals without using 
        any recursion.  Instead we only use a stack of pairs of truth values:

        * Start by pushing True to the stack
        * normal lines go through iff top of Stack is True
        * An IF line is evaluated and truth value pushed on stack
        * An ENDIF pop's the truth value off the stack
        * An ELSE negates top of stack

        All other variables are to produce good error messages.

        Notes
        -----
        To handle ELIF statements we need to push a pair everytime. --- GMK

        """
        stack = [(True, True)]
        iflines = []
        ans = []
        nontrivial = 0
        self.usevars.update(self._init_usevars)

        for self.lineno, line in enumerate(lines, 1):
            state = stack[-1][1]  # current truth value part of top of stack
            m = self._pp_keyword(line)

            # common case first (regular line)
            if (not m):
                # Replace setupvariables with values
                if state:
                    ans.append((self.lineno, line))
                continue

            # Got a PPCmd to handle (not conditional)
            if m in self.normalkeywords:
                # not to be ignored
                if state:
                    self._handle_pp_cmd(m, line)
                continue

            # Now m is one of the conditional commands
            nontrivial = 1
            if m == "IF":
                # need to evaluate an if
                rv = self._eval_cond(line)   # exception passed through to caller
                stack.append((rv, rv))       # push truth value pair to stack
                iflines.append(self.lineno)  # note location of IF
            elif m == "ELSE":
                # switch the most recent state
                a, b = stack[-1]
                b = not a
                stack[-1] = (a, b)
            elif m == 'ELIF':
                # need to evaluate condition
                a, b = stack[-1]
                if not a:
                    # if not already found a true block
                    rv = self._eval_cond(line)  # exception passed through to caller
                    a, b = rv, rv
                else:
                    # ignore this block
                    a, b = True, False
                stack[-1] = a, b
            elif m == "ENDIF":
                # pop the most recent state
                stack.pop()
                iflines.pop()
                if not stack:
                    # stack has become empty
                    msg = ["Configuration Error: ENDIF found without a corresponding IF"]
                    endif_msg = "Extra ENDIF found at line {0} in {1}"
                    endif_msg = endif_msg.format(self.lineno, self.filename)
                    msg.append(endif_msg)
                    msg.append("The real error could be located before the specified line also")
                    raise SetupError("\n".join(msg))
            elif m == "USESETUPVARS":
                # found a use line
                self._check_use(line)
            else:
                # regular line
                raise SetupError("Programming Error")

        # finished processing all lines
        # we have an unbalanced stack
        if len(stack) != 1:
            msg = ("Configuration Error: IF found without a corresponding ENDIF\n"
                   "IF at line {0} in {1} does not have a matching ENDIF")
            msg.format(iflines[-1], self.filename)
            raise SetupError(msg)

        if nontrivial:
            msg = "Processed Version of {0}\n".format(self.filename)
            for lno, line in ans:
                msg += line.replace("\n", "") + '\n'
            warn(msg, RuntimeWarning)

        return ans

    #
    # Special Methods
    #

    def __cmp__(self, other):
        """Alphabetical comparison on unit names (like 'source/io/amr')"""
        if isinstance(other, basestring):
            return cmp(self.unitname, other)
        else:
            return cmp(self.unitname, other.unitname)

    def __repr__(self):
        # Not really a repr string, implement str
        return "Configuration instance for {0}".format(self.unitname)

    #
    # Normal methods follow
    #

    def parent(self):
        """Returns parent unit's full name, empty string if we're orphans"""
        return os.path.split(self.unitname)[0]

    def parse(self):
        """proccess file such that exceptions go all the way to the top."""
        gbls = globals()
        ans = self.process_file()

        for lineno, line in ans:
            rawline = line
            if rawline and rawline[-1] == '\n':
                rawline = rawline[:-1]

            line = strip_comments(line)
            line = line.strip()
            if 0 == len(line):
                continue

            keyword = line.split()[0]
            if keyword not in self.parser_set:
                raise KeyError(self.error_message.format(keyword, self.filename, lineno))

            # Try to match the line with a regex
            self.rawline = rawline  # FIXME this is for STRING PARAM hack
            re_keyword = gbls['RE_' + keyword.upper()]
            mobj = re_keyword.match(line)
            if mobj is not None:
                parse_method_name = "_parse_" + keyword.lower()
                parse_method = getattr(self, parse_method_name)
                parse_method(mobj)
            elif not self.ignore_pp or not (str(msg)[-17:] == " already declared"):
                # Do not abort on multiple occurrence of same PARAMETER in the same file
                # if we are only parsing for documentation purposes (ignore_pp is true) - KW
                m = 'Bad configuation syntax: file {0}, line {1}:\n{2}\n\n{3}'
                m = m.format(self.filename, lineno, rawline, str(msg))
                raise SyntaxError(m)

    def _parse_childorder(self, mobj):
        """CHILDORDER is the order in which some of the children must be processed
        remaining children appear later in lexicographic order"""
        self.childorder = mobj.group(1).split()

    def _parse_datafiles(self, mobj):
        """DATAFILES -> list of wildcard patterns (with absolute path)"""
        self.datafiles.extend([os.path.join(self.unitname, x) for x in mobj.group(1).split()])

    def _parse_d(self, mobj):
        key, comment = mobj.groups()
        if key == '&':
            if self.dkey in self.d:
                self.d[self.dkey] = "{0} {1}".format(self.d[self.dkey], comment)
            else:
                raise SyntaxError('Improper usage of comment continuation')
        else:
            self.dkey = key
            self.d[key] = comment

    def _parse_facevar(self, mobj):
        facevar = mobj.group("varname")
        eosmap = mobj.group("eosmap")

        if not eosmap:
            eosmap = "NONEXISTENT"

        if facevar in self.facevar:
            raise KeyError('FACEVAR {0} already declared'.format(facevar))

        self.facevar[facevar] = (eosmap.upper(), eosmap.upper())

    def _parse_variable(self, mobj):
        variable = mobj.group("varname")
        vartype = mobj.group("type")
        eosmap = mobj.group("eosmap")
        eosmapin = mobj.group("eosmapin")
        eosmapout = mobj.group("eosmapout")
        if eosmapin and eosmapout and eosmap and (eosmap != eosmapin) and (eosmapout != eosmap):
            msg = 'VARIABLE {0} has EOSMAP as well as EOSMAPIN and EOSMAPOUT mappings.'
            msg = msg.format(variable)
            raise ValueError(msg)

        if not eosmapin:
            eosmapin = eosmap

        if not eosmapout:
            eosmapout = eosmap

        if not eosmapin:
            eosmapin = "NONEXISTENT"

        if not eosmapout:
            eosmapout = "NONEXISTENT"

        if not vartype:
            vartype = "GENERIC"

        if variable in self.variable:
            raise KeyError('VARIABLE {0} already declared'.format(variable))

        self.variable[variable] = (vartype.upper(), eosmapin.upper(), eosmapout.upper())

    def _parse_library(self, mobj):
        libname = mobj.group(1).lower()
        libargs = " ".join(mobj.group(2).split())  # trims and removes multiple spaces
        self.library[libname] = libargs

    def _parse_guardcells(self, mobj):
        num = int(mobj.group(1))
        if self.guardcells != None:
            raise AttributeError('GUARDCELLS already declared')
        self.gaurdcells = num

    def _parse_flux(self, mobj):
        self.flux[mobj.group(1)] = 1

    def _parse_scratchvar(self, mobj):
        variable = mobj.group("varname")
        eosmap = mobj.group("eosmap")
        if not eosmap:
            eosmap = "NONEXISTENT"
        self.scratchvar[variable] = (eosmap.upper(), eosmap.upper())

    def _parse_scratchcentervar(self, mobj):
        variable = mobj.group("varname")
        eosmap = mobj.group("eosmap")
        if not eosmap:
            eosmap = "NONEXISTENT"
        self.scratchcentervar[variable] = (eosmap.upper(), eosmap.upper())

    def _parse_scratchfacexvar(self, mobj):
        variable = mobj.group("varname")
        eosmap = mobj.group("eosmap")
        if not eosmap:
            eosmap = "NONEXISTENT"
        self.scratchfacexvar[variable] = (eosmap.upper(), eosmap.upper())

    def _parse_scratchfaceyvar(self, mobj):
        variable = mobj.group("varname")
        eosmap = mobj.group("eosmap")
        if not eosmap:
            eosmap = "NONEXISTENT"
        self.scratchfaceyvar[variable] = (eosmap.upper(), eosmap.upper())

    def _parse_scratchfacezvar(self, mobj):
        variable = mobj.group("varname")
        eosmap = mobj.group("eosmap")
        if not eosmap:
            eosmap = "NONEXISTENT"
        self.scratchfacezvar[variable] = (eosmap.upper(), eosmap.upper())

    def _parse_species(self, mobj):
        name = mobj.group("name")
        n_elec_str = mobj.group("nElec")
        n_elec_str = "" if n_elec_str is None else n_elec_str

        #If the number of electrons are specified in the Config file then
        #nStrElec is not none.  We store the integer representation of this
        #value in the dictionary, or a 1 if this value is not specified.
        if 0 < len(n_elec_str):
            n_elec = int(n_elec_str)
            if n_elec < 1:
                raise ValueError('Must be at least one electron in each species.')
            self.species[name] = n_elec
        else:
            self.species[name] = 1

    def _parse_mass_scalar(self, mobj):
        msname = mobj.group("msname")
        gpname = mobj.group("gpname")
        if gpname:
            gpname = str(gpname).upper()
        eosmap = mobj.group("eosmap")
        eosmapin = mobj.group("eosmapin")
        eosmapout = mobj.group("eosmapout")
        if eosmapin and eosmapout and eosmap and (eosmap != eosmapin) and (eosmapout != eosmap):
            msg = 'MASS_SCALAR {0} has EOSMAP as well as EOSMAPIN and EOSMAPOUT mappings.'
            msg = msg.format(msname)
            raise ValueError(msg)

        if not eosmapin:
            eosmapin = eosmap

        if not eosmapout:
            eosmapout = eosmap

        if not eosmapin:
            eosmapin = "NONEXISTENT"
        if not eosmapout:
            eosmapout = "NONEXISTENT"

        self.mass_scalar[msname] = (gpname, eosmapin.upper(), eosmapout.upper())

    def _parse_particleprop(self, mobj):
        name = mobj.group("name")
        prop_type = mobj.group("type")

        if prop_type == "INTEGER":
            msg = ('Particle properties of type INTEGER are currently not supported, '
                   'use REAL instead to store integer {0} information!')
            msg = msg.format(name)
            raise ValueError(msg)

        if name in self.particleprop:
            raise KeyError('PARTICLEPROP {0} already declared'.format(name))

        self.particleprop[name] = prop_type

    def _parse_particletype(self, mobj):
        """We will store all of the particle type information as a dictionary of tuples.
        The particle type will be key, and the values will be the map & init methods."""
        particle_type = mobj.group("particleType")
        init_method = mobj.group("initMethod")
        map_method = mobj.group("mapMethod")
        adv_method = mobj.group("advMethod")
        if particle_type in self.particletype:
            raise KeyError('PARTICLETYPE {0} already declared'.format(particle_type))
        self.particletype[particle_type] = (init_method, map_method, adv_method)

    def _parse_particlemap(self, mobj):
        name = mobj.group("name")
        type = mobj.group("type")
        varname = mobj.group("varname")
        if name in self.particlemap:
            raise KeyError('PARTICLEMAP {0} already set'.format(name))
        self.particlemap[name] = (type, varname)

    def _parse_suggest(self, mobj):
        sugset = []
        for u in mobj.groups()[0].split():
            if u:
                if u[0] == ".":
                    # using relative path
                    sugset.append(os.path.normpath(os.path.join(self.unitname, u)))
                else:
                    sugset.append(os.path.normpath(u))
        self.suggest.append(sugset)

    def _parse_requires(self, mobj):
        def normpath(path):
            if path.startswith("."):
                path = os.path.join(self.unitname, path)
            return os.path.normpath(path)

        #import pdb; pdb.set_trace()
        # Make a list of all units that followed the "REQUIRES" as delimited by "OR"
        # If only one unit was given, (i.e. there were no "OR"s), re.split() will
        # return a list of one.
        units = re.split("\s+OR\s+", mobj.group("orList"))
        # change any relative paths to absolute paths; remove any empties
        units = [normpath(unit) for unit in units if len(unit) > 0]
        # Note that self.requires is a list of lists
        self.requires.append(units)

    # same as REQUIRES except can be negated by a -without-unit in cmd line
    # also only one unit at a time
    def _parse_requests(self, mobj):
        unit = ""
        u = mobj.group("name")
        if u:
            if u[0] == ".":  # using relative path
                unit = os.path.normpath(os.path.join(self.unitname, u))
            else:
                unit = os.path.normpath(u)

        # do we need to REQUIRE it or not?
        ignore = False  # do not ignore this unit
        reason = "Unknown"
        for x in self.without_units.keys():
            if unit == x:
                ignore = True
                reason = x
            elif unit.startswith(x + os.sep):
                ignore = True
                reason = x

        if not ignore:
            # pretend this was a REQUIRES
            self.requires.append([unit])
            warn("Honoring request for {0}".format(unit), RuntimeWarning)
        else:
            msg = "Ignoring request for {0} (reason: -without-unit={1})".format(unit, reason)
            warn(msg, RuntimeWarning)

    def _parse_parameter(self, mobj):
        name, ptype, pconstant, pvalue, prange = mobj.groups()
        if name in self.parameter:
            raise KeyError('{0} already declared'.format(name))

        if ptype == 'STRING':
            # This is a hack. I "forgot" to accomodate parsers that would need
            # rawlines when I wrote the general framework. 
            pname, pvalue, prange = RE_PSTRING.match(self.rawline).groups()
        else:
            pvalue = pvalue.rstrip()

        self.parameter[name] = (ptype, pvalue, pconstant, prange)

    def _parse_default(self, mobj):
        #self.default = os.path.join(os.path.split(self.filename)[0], 
        self.default = os.path.join(self.unitname, 
                                    os.path.normpath(mobj.group(1)))

    def _parse_exclusive(self, mobj):
        units = mobj.group(1).split()

        # look explicitly for the "*" which means only one directory
        # below can be chosen, ie all are exclusive
        if units[0] == "*":
            name_abspath = os.path.abspath(os.path.split(self.filename)[0])
            dirs = os.listdir(name_abspath)
            for d in dirs:
                if (d[0] != '.') and (os.path.isdir(os.path.join(name_abspath, d))):
                    #units.append(os.path.join(name_abspath, d))
                    units.append(os.path.join(self.unitname, d))
            units.remove("*")
        else:
            #units = [os.path.join(name_abspath, os.path.normpath(x)) for x in units]
            units = [os.path.join(self.unitname, os.path.normpath(x)) for x in units]

        self.exclusive.append(units)

    def _parse_kernel(self, mobj):
        msg = "parse kernel: found {0} in {1}".format(mobj.group(2), self)
        warn(msg, RuntimeWarning)

        if not mobj.group(2):
            self.kernel = self.unitname
        else:
            kname = os.path.join(self.unitname, mobj.group(2))
            self.kernel = kname

        msg = "parse kernel: set to {0}.".format(self.kernel)
        warn(msg, RuntimeWarning)

    def _parse_linkif(self, mobj):
        """LINKIF spaces non-spaces spaces non-spaces spaces ENDOFLINE
        first set of non-spaces is called "filename" and the second set "unitname"."""
        parts = mobj.group(1).split()

        # We are appending a pair (filename, unitname)
        self.linkif.append((os.path.join(self.unitname, parts[0]), parts[1]))

    def _parse_conflicts(self, mobj):
        """CONFLICTS space MODNAME space MODNAME ....
        raise an error if any of MODNAME are also included"""
        for unit in mobj.group(1).split():
            self.conflicts[unit] = 1

    def _parse_ppdefine(self, mobj):
        self.ppdefine[mobj.group("sym")] = mobj.group("val")

    def _parse_nonrep(self, mobj):
        tp = mobj.group("type")
        name = mobj.group("name")
        nlocs = int(mobj.group("nloc"))
        pglob = mobj.group("pglob")
        namef = (mobj.group("namef") or "{name}???").format(name=name)
        
        gbls = globals()
        tpre = gbls['RE_' + tp.upper()]

        parse_method_name = "_parse_" + tp.lower()
        tpparse = getattr(self, parse_method_name)

        locf = 'NONREP_%s_%X_%%06d' % (name, name.__hash__() & 0xffffffff)
        for i in range(nlocs):
            tpparse(tpre.match('{0} {1}'.format(tp, locf % (i + 1))))
        self.nonrep[name] = {'tp': tp,
                             'rpcount': pglob,
                             'nlocs': nlocs,
                             'locf': locf,
                             'namef': namef}

#
# Configuration Union
#

def _update_var_prop(dictout, dictin, infotext):
    """Updates a configuration keyword if that keyword is stored as a dictionary."""
    def check_significant_and_different(l):
        if not l[1] or l[1] == 'NONEXISTENT' or l[1] == 'GENERIC':
            return False
        if not l[0] or l[0] == 'NONEXISTENT' or l[0] == 'GENERIC':
            return False
        return (l[0] != l[1])

    def pick_significant(l):
        if not l[1] or l[1]=='NONEXISTENT' or l[1]=='GENERIC':
            if l[0] and l[0]!='NONEXISTENT' and l[0]!='GENERIC':
                return l[0]
        return l[1]                    
                    

    for k, v in dictin.items():
        if (k not in dictout) or not isinstance(v, collections.Sequence):
            dictout[k] = v
            continue

        # Now handle tuple and list cases
        if not (filter(lambda x: x == 'NONEXISTENT' or x == 'GENERIC', v)):
            zipped_vals = zip(dictout[k], v)
            sig_and_diff = map(check_significant_and_different, zipped_vals)
            cnt_problems = sig_and_diff.count(True)
            first_problem = sig_and_diff.index(True) if 0 < cnt_problems else -1

            # Put some effort into generating meaningful messages.
            if (first_problem >= 0):
                if infotext == "VARIABLE":
                    if first_problem == 0:
                        attr = "TYPE"
                    elif first_problem == 1:
                        if zipped_vals[1][0] == zipped_vals[2][0] and \
                           zipped_vals[1][1] == zipped_vals[2][1]:
                            attr = "EOSMAP"
                        else:
                            attr = "EOSMAPIN"
                    elif first_problem == 2:
                        attr = "EOSMAPOUT"
                    else:
                        attr = "EOSMAP"
                else:
                    if first_problem == 0:
                        if zipped_vals[0][0] == zipped_vals[1][0] and \
                           zipped_vals[0][1] == zipped_vals[1][1]:
                            attr = "EOSMAP"
                        else:
                            attr = "EOSMAPIN"
                    elif first_problem==1:
                        attr = "EOSMAPOUT"
                    else:
                        attr = "EOSMAP"

                msg = 'Conflicting specifications for {0} {1} {2}: "{3}" and "{4}".'
                msg = msg.format(infotext, k, attr, zipped_vals[first_problem][0], 
                                                    zipped_vals[first_problem][1])
                raise ValueError(msg)
            dictout[k] = v
        else:
            zipped_vals = zip(dictout[k], v)
            sig_and_diff = map(check_significant_and_different, zipped_vals)
            cnt_problems = sig_and_diff.count(True)
            first_problem = sig_and_diff.index(True) if 0 < cnt_problems else -1

            # Put some effort into generating meaningful messages.
            if (first_problem >= 0):
                if infotext == "VARIABLE":
                    if first_problem == 0:
                        attr = "TYPE"
                    elif first_problem == 1:
                        if zipped_vals[1][0] == zipped_vals[2][0] and \
                           zipped_vals[1][1] == zipped_vals[2][1]:
                            attr = "EOSMAP"
                        else:
                            attr = "EOSMAPIN"
                    elif first_problem == 2:
                            attr = "EOSMAPOUT"
                    else:
                            attr = "EOSMAP"
                else:
                    if first_problem==0:
                        if zipped_vals[0][0] == zipped_vals[1][0] and \
                           zipped_vals[0][1] == zipped_vals[1][1]:
                            attr = "EOSMAP"
                        else:
                            attr = "EOSMAPIN"
                    elif first_problem==1:
                        attr = "EOSMAPOUT"
                    else:
                        attr = "EOSMAP"
                msg = 'Conflicting specifications for %s %s %s: "%s" and "%s".'
                msg = msg.format(infotext, k, attr, zipped_vals[first_problem][0],
                                                    zipped_vals[first_problem][1])
                raise ValueError(msg)
            merged_vals = map(pick_significant, zipped_vals)
            dictout[k] = merged_vals
                    

class ConfigurationUnion(object):
    """Container for collation of Configuration instances."""

    def __init__(self, configs):
        """Simply takes a list of config files.

        Parameters
        ----------
        configs : sequence
            Configuration instances.

        """
        self.datafiles = []
        self.facevar = {}
        self.flux = {}
        self.scratchvar = {}
        self.scratchcentervar = {}
        self.scratchfacexvar = {}
        self.scratchfaceyvar = {}
        self.scratchfacezvar = {}
        self.guardcells = None
        self.library = {}
        self.particleprop = {}
        self.particletype = {}
        self.particlemap = {}
        self.species = {}
        self.mass_scalars = {}
        self.mass_scalar_groups = {}  # list of names of groups (in arbitrary order)
        self.linkif = []
        self.variable = {}
        self.ppdefines = {}
        self.nonrep = {}

        self.particles_unit_included = False

        self.collate(configs)
        self.set_defaults()
        self.simplify()

    #
    # Normal methods follow
    #

    def collate(self, configs):
        """Collates info from an ordered list of Configuration instances.

        Parameters
        ----------
        configs : sequence
            Configuration instances.

        """
        partproplocs = {}  # for err msg
        speclocs = {}      # dictionary storing which units define which species.
        mscalarlocs = {}   # for err msg

        for config in configs:
            # check if we have included any particle units.
            # We only check for Particles subdirectories because the top level
            # Particles directory (full of stubs) is included in ALL simulations.
            if config.unitname.startswith("Particles/"):
                self.particles_unit_included = True

            # self.library is a dictionary of lists of arguments
            # In perfect world all entries of this list will be the same
            for lib in config.library.keys():
                if lib not in self.library:
                    self.library[lib] = []
                self.library[lib].append(config.library[lib])

            _update_var_prop(self.scratchvar, config.scratchvar, "SCRATCHVAR")
            _update_var_prop(self.scratchcentervar, config.scratchcentervar, "SCRATCHCENTERVAR")
            _update_var_prop(self.scratchfacexvar, config.scratchfacexvar, "SCRATCHFACEXVAR")
            _update_var_prop(self.scratchfaceyvar, config.scratchfaceyvar, "SCRATCHFACEYVAR")
            _update_var_prop(self.scratchfacezvar, config.scratchfacezvar, "SCRATCHFACEZVAR")
            _update_var_prop(self.facevar, config.facevar, "FACEVAR")
            self.flux.update(config.flux)
            _update_var_prop(self.variable, config.variable, "VARIABLE")
            self.ppdefines.update(config.ppdefine)
            self.datafiles.extend(config.datafiles)

            self.nonrep.update(config.nonrep)

            for name, (group, eos1, eos2) in config.mass_scalar.items():
                if name in self.mass_scalars:
                    currval = self.mass_scalars[name]
                    gotval = True
                else:
                    self.mass_scalars[name] = (group, eos1, eos2)
                    gotval = False

                # different groups specified
                if gotval and (currval != group):
                    msg = ['MASS SCALAR [{0}] has two groups '.format(name)]
                    msg.append('    {0} in {1}'.format(currval, config.unitname))
                    msg.append('    {0} in {1}'.format(group, mscalarlocs[name]))
                    raise ValueError("\n".join(msg))

                # everything is fine now
                mscalarlocs[name] = config.unitname  # store where it was defined
                if group:
                    self.mass_scalar_groups[group] = 1

            # Check whether we have already encountered a species of the same
            # name.  If so, let the user know which Config files contain
            # the same SPECIES definition.  Otherwise, add the species to
            # the main dictionary.
            for (key_string, num_elec) in config.species.items():
                if key_string in self.species:
                    # We have defined the same species in multiple Config files.
                    msg = ['SPECIES {0} defined twice'.format(key_string)]
                    msg.append('    {0} defines the species {1}'.format(config.name, config.species))
                    msg.append('    {0} defines the species {1}'.format(speclocs[key_string],
                                                                      self.species))
                    raise ValueError("\n".join(msg))
                else:
                    self.species[key_string] = num_elec     # Copy species into main dictionary.
                    speclocs[key_string] = config.unitname  # Store config name that defines species.

            for ud in config.linkif:
                if ud not in self.linkif:
                    self.linkif.append(ud)

            if config.guardcells != None:
                if self.guardcells == None:
                    self.guardcells = config.guardcells
                else:
                    self.guardcells = max(config.guardcells, self.guardcells)

            for (prop, prop_type) in config.particleprop.items():
                if config.particleprop.get(prop, prop_type) != prop_type:
                    msg = ['Integer/Real PARTICLEPROP MISMATCH with %s '.format(prop)]
                    msg.append('    has type %s in %s'.format(prop_type, config.name))
                    msg.append('    has type %s in %s'.format(self.particleprop[prop],
                                                              partproplocs[prop]))
                    raise SetupError("\n".join(msg))
                else:
                    self.particleprop[prop] = prop_type
                    partproplocs[prop] = config.unitname

            if config.particletype != None:
                for particle_type in config.particletype.keys():
                    if particle_type in self.particletype:
                        raise SetupError('PARTICLETYPE {0} already declared'.format(particle_type))
                #We've ensured we have no common keys, so copy the data across:
                self.particletype.update(config.particletype)

            for prop, map in config.particlemap.items():
                if config.particlemap.get(prop, map) != map:
                    msg = ['PARTICLEMAP  mismatch with {0} '.format(prop)]
                    raise SetupError("\n".join(msg))
                else:
                    self.particlemap[prop] = map

    def set_defaults(self):
        if self.guardcells is None:
            self.guardcells = 4
            warn('number of guard cells not specified, using default GUARDCELLS 4', 
                 RuntimeWarning)

    def simplify(self):
        """Store stuff that can be computed from configuration info.
        We'll use lowercase keys to separate things.

        In many cases self['UPPERCASE'] is a dictionary. We define
        self['uppercase'] as its sorted list of keys. """

        ppd = self.ppdefines  # Speeds up lookup time
        ppds = []
        for x in ppd.keys():
            if ppd[x]:
                ppds.append("#define {0} {1}".format(x.upper(), ppd[x]))
            else:
                ppds.append("#define {0}".format(x.upper()))

        ppds.sort()
        self.ppdefines_lines = ppds

        # need this since order in which a dict returns it's keys is not determined.
        self.variable_names = self.variable.keys()
        self.variable_names.sort()
        tmp_list = [self.variable[name] for name in self.variable_names]
        self.var_types = [x for (x, y, z) in tmp_list]          # list of corresponding TYPES
        self.eosmapin_unkvars = [y for (x, y, z) in tmp_list]   # list of Eos_map roles
        self.eosmapout_unkvars = [z for (x, y, z) in tmp_list]  # list of Eos_map roles
        self.nparticletypes = len(self.particletype)

        # Add a quick sense check.  Check that a simulation having 0 particle
        # types really does not use particles.
        particle_types_included = (self.nparticletypes > 0)
        if (particle_types_included != self.particles_unit_included):
            msg = ('ERROR! Particle conflict\nHave we included a particle config: {0}.\n'
                   'Have we included a particle type: {0}.')
            msg = msg.format(self.particles_unit_included, particle_types_included)
            raise ValueError(msg)

        # When we define multiple particle types we must add a field named
        # TYPE_PART_PROP so FLASH can distinguish between particle types.
        if 'type' in self.particleprop:
            raise SyntaxError('ERROR! TYPE_PART_PROP is a reserved definition\n')
        else:
            if self.nparticletypes > 1:
                self.particleprop['type'] = 'REAL'

        self.intproperty = []
        self.realproperty = []
        pkeys = self.particleprop.keys()

        # DEV: not using intproperty right now
        self.intproperty = [x for x in pkeys if self.particleprop[x] == 'INTEGER']
        self.intproperty.sort()
        if self.intproperty:
            self.n_int_props = len(self.intproperty)
        else:
            self.n_int_props = 1

        self.realproperty = [x for x in pkeys if self.particleprop[x] == 'REAL']
        self.realproperty.sort()
        if self.realproperty:
            self.n_real_props = len(self.realproperty)
        else:
            self.n_real_props = 1

        # Provide a custom compare function to the list's sort method.
        # Sorts list by the number of electrons (field 1) and then alphabetically (field 0).
        # The second condition is only evaluated if the first condition is false.
        list_tuples = self.species.items()
        list_tuples.sort(lambda x, y: cmp(int(x[1]), int(y[1])) or cmp(x[0], y[0]))

        total_species = 0
        self.species_names = []
        for (element, num_electrons) in list_tuples:

            #We must sum over all electrons in all species.
            if num_electrons > 1:
                total_species = total_species + (num_electrons + 1)
                len_appended_number = len(str(num_electrons))
            else:
                total_species += 1
                len_appended_number = 0

            # There is an inbuilt 4-character limit in Paramesh which truncates names.
            # Warn the user if they have exceeded this limit.
            # If there is only one ion for a particular element then we append no number.
            if (len(element) + len_appended_number) > 4:
                raise SetupError('Constructed element name for "{0}" will '
                                 'exceed 4 character limit'.format(element))

            # Place the unique name in the list.
            self.species_names.append(element)
            if num_electrons > 1:
                for electron in range(1, num_electrons + 1, 1):
                    self.species_names.append("{0}{1}".format(element, electron))

        self.nspecies = total_species
        self.nmassscalars = len(self.mass_scalars)
        self.massscalars = self.mass_scalars.keys()
        self.massscalars.sort()
        tmp_list = [self.mass_scalars[var] for var in self.massscalars]

        self.eosmapin_ms = [y for (x, y, z) in tmp_list]
        self.eosmapout_ms = [z for (x, y, z) in tmp_list]

        msg = self.mass_scalar_groups.keys()
        msg.sort()
        self.massscalars_map = {}
        self.massscalars_group_map = {}

        for name, (group, eos1, eos2) in self.mass_scalars.items():
            if not group:
                continue

            # index return 0 based, our groups start at 1
            self.massscalars_map[name] = msg.index(group) + 1

            # index return 0 based, our groups start at 1
            self.massscalars_group_map[group] = msg.index(group) + 1

        # now masscalars_map is a dictionary mapping scalar to its group number
        # only mass_scalars which need to be renormed appear in this map
        # mass_scalars_group_map maps group names to numbers

        self.nfacevars = len(self.facevar)
        self.facevar_names = self.facevar.keys()
        self.facevar_names.sort()
        #self.eos_facevars = [ self.facevar[var] for var in self.facevar ] # list of Eos_map roles.
        tmp_list = [self.facevar[var] for var in self.facevar_names]
        self.eosmapin_facevars = [x for (x, y) in tmp_list]   # list of Eos_map roles.
        self.eosmapout_facevars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.flux_names = self.flux.keys()
        self.flux_names.sort()

        self.scratchvar_names = self.scratchvar.keys()
        self.scratchvar_names.sort()
        #self['eos_scratchvars'] = [ self.scratchvar[var] for var in self.scratchvar ] # list of Eos_map roles.
        tmp_list = [self.scratchvar[var] for var in self.scratchvar_names]
        self.eosmapin_scratchvars = [x for (x, y) in tmp_list]   # list of Eos_map roles.
        self.eosmapout_scratchvars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.scratchcentervar_names = self.scratchcentervar.keys()
        self.scratchcentervar_names.sort()
        tmp_list = [self.scratchcentervar[var] for var in self.scratchcentervar_names]
        self.eosmapin_scratchcentervars = [x for (x, y) in tmp_list]   # list of Eos_map roles.
        self.eosmapout_scratchcentervars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.scratchfacexvar_names = self.scratchfacexvar.keys()
        self.scratchfacexvar_names.sort()
        tmp_list = [self.scratchfacexvar[var] for var in self.scratchfacexvar_names]
        self.eosmapin_scratchfacexvars = [x for (x, y) in tmp_list]   # list of Eos_map roles.
        self.eosmapout_scratchfacexvars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.scratchfaceyvar_names = self.scratchfaceyvar.keys()
        self.scratchfaceyvar_names.sort()
        tmp_list = [self.scratchfaceyvar[var] for var in self.scratchfaceyvar_names]
        self.eosmapin_scratchfaceyvars = [x for (x, y) in tmp_list]   # list of Eos_map roles.
        self.eosmapout_scratchfaceyvars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.scratchfacezvar_names = self.scratchfacezvar.keys()
        self.scratchfacezvar_names.sort()
        tmp_list = [self.scratchfacezvar[var] for var in self.scratchfacezvar_names]
        self.eosmapin_scratchfacezvars = [x for (x, y) in tmp_list]  # list of Eos_map roles.
        self.eosmapout_scratchfacezvars = [y for (x, y) in tmp_list]  # list of Eos_map roles.

        self.nvar = self.nspecies + self.nmassscalars + len(self.variable)
        self.nflux = self.nspecies + self.nmassscalars + len(self.flux)
        self.max_plot_vars = self.nvar  # + (self.nfacevars*3) + self.nflux

        self.particle_type = []
        self.init_method = []
        self.mapmethod = []
        self.advmethod = []

        # We create a sorted temporary list "tmp_list", which we then copy into
        # a new list "particle_list" with passive particles at the start.
        tmp_list = self.particletype.items()
        tmp_list.sort(lambda x, y: cmp(x[0], y[0]))  # Sort by particle type.

        i = 0
        particle_list = []
        for (particle_type, (init_method, map_method, adv_method)) in tmp_list:
            if (particle_type.upper() == "PASSIVE"):
                particle_list.append(tmp_list[i])
                del tmp_list[i]
                break
            i = i + 1
        particle_list.extend(tmp_list)

        for (particle_type, (init_method, map_method, adv_method)) in particle_list:
            self.particle_type.append(particle_type.upper())
            self.init_method.append(init_method.upper())
            self.mapmethod.append(map_method.upper())
            self.advmethod.append(adv_method.upper())

        # build separate lists for the different grid-variable typed particle maps
        pmkeys = self.particlemap.keys()

        self.particlemaps_variable = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'VARIABLE']
        self.particlemaps_species = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'SPECIES']
        self.particlemaps_mscalar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'MASS_SCALAR']
        self.particlemaps_facex = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'FACEX']
        self.particlemaps_facey = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'FACEY']
        self.particlemaps_facez = [(x, self.particlemap[x][1]) for x in pmkeys \
                                        if self.particlemap[x][0] == 'FACEZ']
        self.particlemaps_scratchvar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                            if self.particlemap[x][0] == 'SCRATCHVAR']
        self.particlemaps_scratchcentervar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                                if self.particlemap[x][0] == 'SCRATCHCENTERVAR']
        self.particlemaps_scratchfacexvar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                                if self.particlemap[x][0] == 'SCRATCHFACEXVAR']
        self.particlemaps_scratchfaceyvar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                                if self.particlemap[x][0] == 'SCRATCHFACEYVAR']
        self.particlemaps_scratchfacezvar = [(x, self.particlemap[x][1]) for x in pmkeys \
                                                if self.particlemap[x][0] == 'SCRATCHFACEZVAR']

        self.nparticlemaps_variable = len(self.particlemaps_variable)
        self.nparticlemaps_facex = len(self.particlemaps_facex)
        self.nparticlemaps_facey = len(self.particlemaps_facey)
        self.nparticlemaps_facez = len(self.particlemaps_facez)
        self.nparticlemaps_scratchvar = len(self.particlemaps_scratchvar)
        self.nparticlemaps_scratchcentervar = len(self.particlemaps_scratchcentervar)
        self.nparticlemaps_scratchfacexvar = len(self.particlemaps_scratchfacexvar)
        self.nparticlemaps_scratchfaceyvar = len(self.particlemaps_scratchfaceyvar)
        self.nparticlemaps_scratchfacezvar = len(self.particlemaps_scratchfacezvar)

        #check the consistency of the particle maps
        for part_var, grid_var in self.particlemaps_variable:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.variable_names.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {0}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_species:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.species.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_mscalar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.massscalars.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_facex:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.facevar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_facey:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.facevar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_facez:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.facevar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_scratchvar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.scratchvar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_scratchcentervar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.scratchcentervar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_scratchfacexvar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.scratchfacexvar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_scratchfaceyvar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.scratchfaceyvar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

        for part_var, grid_var in self.particlemaps_scratchfacezvar:
            if (self.realproperty.count(part_var) != 1):
                msg = 'fatal:  bad particle map for {0}; no such property'
                msg = msg.format(part_var)
                raise ValueError(msg)

            if (self.scratchfacezvar.count(grid_var) == 0):
                msg = 'fatal:  bad particle map for {0}; no such variable {1}'
                msg = msg.format(part_var, grid_var)
                raise ValueError(msg)

