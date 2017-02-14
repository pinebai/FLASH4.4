import os
import re
from warnings import warn

# relative imports needed!  
from ..dsl.configuration import Configuration, ConfigurationUnion
from . import setup_globals
from .setup_globals import gvars, SetupError
from .utils import search_paths
from .rp_info import RPInfo
from .var_info import VarInfo
from .template import Template

#
# Configuration List
#

class ConfigurationList(object):
    """Declares a class ConfigurationtList which contains a list of unit names
    together with methods to manipulate them.  This must be aware of how the 
    problem is being setup.
    """

    def __init__(self, ignore_pp=False):
        self.units = {}  # dictionary mapping name to FlashUnit object
        self.ignore_pp = ignore_pp

    def check_empty(self):
        # programming error
        if self.units:
            gvars.out.put("Programming Error: self.units should be empty here",
                          setup_globals.ERROR)
            raise SetupError

    def remove_units(self, unitnames):
        """Remove these units also
        returns the number of units actually removed"""
        if isinstance(unitnames, basestring):
            unitnames = [unitnames]

        count = 0
        for unitname in unitnames:
            if not unitname:
                continue
            if unitname in self.units:
                count += 1
                del self.units[unitname]

        return count

    def remove_subtrees(self, unitnames):
        """remove these units and all their children"""
        # speed up common case
        if not unitnames:
            return

        rmunitnames = []
        for unitname in unitnames:
            if unitname.endswith(os.sep):
                rmunitnames.append(unitname[:-1])
            else:
                rmunitnames.append(unitname)

        # rmunitnames has all unitnames not ending with "/"
        rmlist = []  # names of units to be removed
        for unitname in rmunitnames:
            usep = unitname + os.sep
            for uname in self.units.keys():
                if uname == unitname or uname.startswith(usep):
                    rmlist.append(uname)

        for uname in rmlist:
            del self.units[uname]

    def add_units(self, unitnames, comment=""):
        """Add these units and returns the number of units actually added."""
        if isinstance(unitnames, basestring):
            unitnames = [unitnames]

        count = 0
        for unitname in unitnames:
            if 0 == len(unitname):
                continue

            if unitname in self.units:
                continue

            source_abspath = [os.path.abspath(sd) for sd in gvars.source_path]
            if os.path.isabs(unitname) and unitname in source_abspath:
                continue

            if count == 0 and comment:
                gvars.out.put(comment, setup_globals.INFO)
                gvars.out.push()

            count += 1
            if comment:
                gvars.out.put(unitname, setup_globals.IMPINFO)

            top_unit = (unitname in gvars.topUnitNames)

            # find and add unit
            name = os.path.normpath(unitname)  # something like 'IO/common/hdf5/'
            name_path = name if name.startswith(gvars.project_simulations_dir) \
                             else search_paths(gvars.source_path, name)
            name_abspath = os.path.abspath(name_path)

            if not os.path.isdir(name_abspath):
                raise SetupError('Unit {0} not found (at {1})'.format(name, name_abspath))

            conf_names = ['config.py', 'Config.py', 'Config']
            confs_present = [f for f in os.listdir(name_abspath) if f in conf_names]
            config_filename = None if 0 == len(confs_present) \
                                   else os.path.join(name_abspath, confs_present[0])
            #if config_filename is None:
            #    continue

            if config_filename is not None and config_filename.endswith('.py'):
                # FIXME do something special to handle python config files
                conf = None
            else:
                conf = Configuration(config_filename, top_unit=top_unit, unitname=name,
                                     usevars=gvars.setup_vars.vars, ignore_pp=self.ignore_pp,
                                     without_units=gvars.withoutUnits, unitpath=name_abspath)

            self.units[unitname] = conf

        if count > 0 and comment:
            gvars.out.pop()

        return count

    def check_top_units(self):
        for top_unit in gvars.topUnitNames:
            if top_unit not in self.units:
                continue

            init_file = os.path.basename(top_unit) + "_init.F90"
            if not os.path.isfile(os.path.join(self.units[top_unit].filename, init_file)):
                msg = ("WARNING: Missing 'init' File\n"
                       "Each API-level unit directory must "
                       "directly contain an initialization file\n"
                       'whose name matches the name of that directory plus "_init.F90"\n'
                       'The directory "{0}" does not contain a file "{1}"')
                msg = msg.format(top_unit, init_file)
                warn(msg, RuntimeWarning)

    def add_top_units(self):
        self.add_units(gvars.topUnitNames)

    def check_suggest(self):
        """list of units which have been suggested but not included."""
        badunits = []
        unames = self.units.keys()
        for uname in unames:
            for suglist in self.units[uname].suggest:
                badlist = 1
                for sug in suglist:
                    if sug in unames:
                        badlist = None

                # append the alternatives
                if badlist:
                    badunits.append((uname, suglist))

        # some suggestion not taken
        if badunits:
            gvars.out.put("\nIGNORED SUGGESTIONS", setup_globals.IMPINFO)
            for (uname, suglist) in badunits:
                gvars.out.put("{0} suggested one of the following units as well: ".format(uname),
                              setup_globals.IMPINFO)
                gvars.out.push()
                for sug in suglist:
                    gvars.out.put(sug, setup_globals.IMPINFO)
                gvars.out.pop()
            gvars.out.put("", setup_globals.IMPINFO)

    def check_requirements(self):
        for unit in self.units.values():
            # will be a list of lists
            # each member of this list-within-a-list will represent
            # a series of alternative Units as specified in a Config
            # file by the syntax "REQUIRES A OR B OR C".
            sets_of_alternatives = unit.requires
            for set_of_alternatives in sets_of_alternatives:
                for unit_name in set_of_alternatives:
                    if unit_name in self.units:
                        # We have one in the list of required units
                        # so don't bother checking the others
                        break
                else:
                    raise SetupError('{0} REQUIRES {1}, not included'.format(unit.name,
                                     " or ".join(set_of_alternatives)))

    def check_exclusivity(self):
        gvars.out.put("Checking if Exclusive units are included", setup_globals.DEBUG)

        # units contains everybody's parents
        for unit in self.units.values():
            for group in unit.exclusive:
                # No two elements of group must be in units
                a = None
                for b in group:
                    if b not in self.units:
                        continue
                    if a:
                        raise SetupError('{0} and {1} are exclusive'.format(a, b))
                    a = b

    def check_conflicts(self):
        gvars.out.put("Checking for Conflicting units", setup_globals.DEBUG)
        for unit in self.units.values():
            for b in unit.conflicts:
                if b in self.units:
                    raise SetupError("setup Error: requested unit "
                                     "{0} CONFLICTS with {1}".format(b, unit))

    def get_unit_names_from_file(self, file):
        if not os.path.isfile(file):
            gvars.out.put('cannot access {0} file'.format(file), setup_globals.ERROR)
            gvars.out.push()
            gvars.out.put('* Either use the -auto option or', setup_globals.ERROR)
            gvars.out.put('* specify the Units file using the -unitsfile=<filename> option',
                          setup_globals.ERROR)
            gvars.out.pop()
            raise SetupError('No Units file found', "NOUNITS")

        gvars.out.put('scanning {0} file for included units'.format(file), setup_globals.IMPINFO)
        ans = []
        for line in open(file).readlines():
            rawline = line
            if line.count('#'):
                line = line[:line.find('#')]
            line = line.strip()

            if not line:
                continue

            try:
                a, b = line.split()
                if a != 'INCLUDE':
                    raise SetupError
            except (ValueError, SetupError):
                raise SetupError('Bad syntax:\n' + rawline)

            ans.append(b)

        return ans

    def add_units_from_file(self):
        fname = os.path.join(gvars.project_setup_dir, setup_globals.UNITS_FILENAME)
        self.add_units(self.get_unit_names_from_file(fname), "Scanning Units file")

    def get_default_units(self):
        """return names of default units of units already present"""
        return [unit.default for unit in self.units.values()]

    def add_default_units(self):
        """add default units of all units"""
        while self.add_units(self.get_default_units(), "") > 0:
            pass

    def add_kernel_units(self):
        """check to see which units are Kernel units
        if we have a Kernel unit then we want to
        include all units and directories under the
        unit specified to be a Kernel
        """
        ans = []
        gvars.out.push()
        gvars.out.put('*** KERNEL *** add_kernel_units *** START ***', setup_globals.DEBUG)
        for unit in self.units.values():
            if unit.kernel:
                gvars.out.put('Processing KERNEL for {0}'.format(unit), setup_globals.DEBUG)
                gvars.out.push()
                kdir = search_paths(gvars.source_path, unit.kernel)
                if (kdir is None) or not os.path.isdir(kdir):
                    msg = "Invalid KERNEL {0}, not a directory, requested by {1}.\n"
                    msg = msg.format(unit.kernel, unit)
                    raise SetupError(msg)
                gvars.out.pop()
                # FIXME recursive_get_dir() needs to be run from the source dir...
                tmpdir = kdir[:kdir.index(unit.kernel)]
                cwd = os.getcwd()
                os.chdir(tmpdir)
                ans.extend(self.recursive_get_dir(unit.kernel))
                os.chdir(cwd)
        self.add_units(ans)
        gvars.out.put('*** KERNEL *** add_kernel_units *** END ***', setup_globals.DEBUG)
        gvars.out.pop()

    def add_grid_interpolation(self):
        """Different interpolation units need to be added depending on which
        Grid has been chosen. This is a somewhat hacky solution to that
        problem. As far as I know this is the only instance of a command-
        line option whose effect depends on which other Units have already
        been included. -nttaylor
        """
        # Note that no interpolation Unit needs to be added for UG
        ans = []
        if gvars.gridInterpolation == setup_globals.GRID_INTERP_MONOTONIC:
            for unit_name in self.units.keys():
                if unit_name.startswith("Grid/GridMain/paramesh/paramesh4"):
                    # Check if added this already
                    for unit_name2 in self.units.keys():
                        if unit_name2.startswith("Grid/GridMain/paramesh/interpolation/Paramesh4"):
                            break
                    else:
                        # add the required interpolation unit
                        ans.append("Grid/GridMain/paramesh/interpolation/Paramesh4")
                    break
                elif unit_name.startswith("Grid/GridMain/paramesh/Paramesh2"):
                    # Check if added this already
                    for unit_name2 in self.units.keys():
                        if unit_name2.startswith("Grid/GridMain/paramesh/Paramesh2/monotonic"):
                            break
                    else:
                        # add the required interpolation unit
                        ans.append("Grid/GridMain/paramesh/Paramesh2/monotonic")
                    break
        elif gvars.gridInterpolation == setup_globals.GRID_INTERP_NATIVE:
            # This is the only other possiblity, and all it means is that we
            # don't have to add any special units, so there's nothing to do.
            pass

        self.add_units(ans)

    def recursive_get_dir(self, path):
        def vfunc(ans, dname, fnames):
            gvars.out.put('...while walking, ans='
                          '{0}, dname="{1}", fnames={1} ...'.format(ans, dname, fnames),
                          setup_globals.DEBUG)
            dont_descend = []

            # remove all (a) ".files" and (b) non-directories from fnames
            for x in fnames:
                if x[0] == ".":
                    # names starting with . to be ignored
                    dont_descend.append(x)
                    continue

                jname = os.path.join(dname, x)
                if not os.path.isdir(jname):
                    # not a directory, also ignored
                    dont_descend.append(x)
                    continue

            # removal in place so recursion does not go there
            for x in dont_descend:
                fnames.remove(x)

            # if we descended into this dname directory, append it to answer list
            ans.append(dname)

        ans = []
        gvars.out.put('Will walk {0} ...'.format(path), setup_globals.DEBUG)
        os.path.walk(path, vfunc, ans)

        # now ans contains the list of all directories inside path
        # except "." directories and their children
        gvars.out.put('...and found {0} by walking.'.format(ans), setup_globals.DEBUG)
        return ans

    def parent_units(self):
        return [unit.parent() for unit in self.units.values()]

    def add_parent_units(self):
        while self.add_units(self.parent_units()) > 0:
            pass

    def remove_parent_units(self):
        self.remove_units(self.parent_units())

    def get_required_units(self):
        """list of unitnames to be added"""
        ans = []
        for unit in self.units.values():
            # will be a list of lists
            # each member of this list-within-a-list will represent
            # a series of alternative Units as specified in a Config
            # file by the syntax "REQUIRES A OR B OR C". Only one of
            # these units need be appended to 'ans'
            sets_of_alternatives = unit.requires
            for set_of_alternatives in sets_of_alternatives:
                for unit_name in set_of_alternatives:
                    if unit_name in self.units:
                        # we already have one in the list of required units
                        # so don't bother checking the others
                        break
                else:
                    # only if we didn't have any of the listed alternatives
                    # do we choose the first one by appending it to 'ans'
                    ans.append(set_of_alternatives[0])
        return ans

    def add_required_units(self):
        """Satisfy requirements"""
        while self.add_units(self.get_required_units()) > 0:
            pass

    def write_units_file(self):
        """Disable old method of object dir creation, write to current dir"""
        outfile = os.path.join(gvars.project_setup_dir, setup_globals.UNITS_FILENAME)
        if os.path.isfile(outfile):
            gvars.out.push()
            gvars.out.put('Backing up {0}'.format(outfile), setup_globals.INFO)
            os.rename(outfile, outfile + '.bak')
            gvars.out.pop()

        outfd = open(outfile, 'w')
        outfd.write('#Units file for {0} generated by setup \n\n'.format(gvars.simulation_name))

        # Dont print parent units (makes file concise)
        p_units = self.parent_units()
        c_units = [x for x in self.units.keys() if x not in p_units]
        c_units.sort()
        for unitname in c_units:
            #if (gvars.project_base_dir is not None and 
            #    unitname.startswith(gvars.project_base_dir)):
            #    unitname = os.path.relpath(unitname, gvars.project_base_dir)
            outfd.write('INCLUDE {0}\n'.format(unitname))
        outfd.close()

    #
    # User calls methods below
    #

    def get_list(self):
        list = self.units.keys()
        list.sort()
        return list

    def get_link_order(self):
        sim_root_dir = setup_globals.SIM_SETUP_DIR.split("/")[0]

        def earlier(corder, childa, childb):
            """Return true if child "unita" should be linked earlier"""
            # code comes here only if we have a non-trivial corder
            if childa in corder:
                try:
                    return cmp(corder.index(childa), corder.index(childb))
                except:
                    # unitb is not in the list
                    return -1
            else:
                # unita not present in list
                if childb in corder:
                    return 1
                else:
                    # both not in corder
                    return cmp(childa, childb)

        def defcmp(childa, childb):
            """default comparison between different children"""
            # Push simulations to the end
            if childa == sim_root_dir:
                return 1
            if childb == sim_root_dir:
                return -1

            # normal compare
            return cmp(childa, childb)

        def apply_aliases(apath):
            ap = "/".join(apath)
            apaths = [ap.replace(k, v) for k, v in gvars.unit_aliases.items() 
                                                if ap.startswith(k)]
            if len(apaths) == 0:
                return apath
            elif len(apaths) == 1:
                return apaths[0].split('/')
            else:
                raise ValueError("Too many matching aliases for {0}: {1}".format(ap, apaths))

        def compare(apath, bpath):
            """return 0 if equal, -1 if patha < pathb, 1 if patha > pathb"""
            if apath == bpath:
                return 0

            # Apply aliases for comparisons
            apath = apply_aliases(apath)
            bpath = apply_aliases(bpath)

            # find youngest common ancestor
            common = []
            patha = apath[:]  # copy is essential here
            pathb = bpath[:]
            while patha and pathb and (patha[0] == pathb[0]):
                common.append(patha[0])
                del patha[0]
                del pathb[0]
            parent = "/".join(common)  # found youngest common parent
            # patha is a prefix of pathb
            if not patha:
                return -1

            # pathb is a prefix of patha
            if not pathb:
                return 1

            # Now need to compare childa and childb according to order specified by parent
            childa = patha[0]
            childb = pathb[0]

            # if trouble accessing parent's order (pick default ordering)
            try:
                corder = self.units[parent].childorder
                if corder:
                    return earlier(corder, childa, childb)
                else:
                    return defcmp(childa, childb)
            except:
                return defcmp(childa, childb)

        # returns names of units present in a fixed order
        # split unitnames into list describing path
        pathlist = [unitname.split("/") for unitname in self.units.keys()]

        # First do a lexicographic sort as that should get most of elements in place
        # since this is optimized in python, we have done most of the work without using
        # our fancy ordering
        pathlist.sort()

        # now use our fancy sorting routine, so hopefully there are few calls to our
        # comparator algorithm
        pathlist.sort(compare)

        # combine paths back to strings
        pl = ["/".join(x) for x in pathlist]  # contract paths to unitnames
        return pl

    def get_config_info(self, **kw):
        """Generate FlashUnitUnion class based on given units"""
        cinfo = ConfigurationUnion(self.units.values(), **kw)
        cinfo.unit_names = self.get_list()
        return cinfo

    def generate_units_file(self):
        """called only when '-auto' flag is passed"""
        gvars.out.put('generating default Units file', setup_globals.IMPINFO)
        self.check_empty()

        # add the appropriate unit from simulation
        self.add_units(gvars.simulation_dir)
        if gvars.simulation_dir.startswith(gvars.project_simulations_dir):
            # if in a local simultion dir, we must manually add the SimulationMain unit
            self.add_units(setup_globals.SIM_SETUP_DIR)

        # add units specified on command-line with "--with-unit"
        self.add_units(gvars.withUnits)

        while True:
            oldnames = self.get_list()  # names of all units in this FlashUnitList instance

            self.add_parent_units()
            self.add_required_units()
            self.remove_parent_units()

            # start over if we added any required units
            if oldnames != self.get_list():
                continue

            self.add_default_units()
            self.add_kernel_units()

            self.add_grid_interpolation()

            # Nothing new added so we are done
            if oldnames == self.get_list():
                break

        # now kill units
        if gvars.killUnits:
            gvars.out.put('Killing subtrees: {0}'.format(
                          ",".join(gvars.killUnits.keys())),
                          setup_globals.INFO)  # USER BEWARE
            self.remove_subtrees(gvars.killUnits.keys())

        # Tell User about the list of units
        gvars.out.push()
        for unitname in self.get_list():
            gvars.out.put(unitname, setup_globals.IMPINFO)
        gvars.out.pop()

        self.write_units_file()

        # If user had not used "-auto", we would not even be in this branch
        # of code, so now that the Units file is written, we clear out the
        # units dictionary again. Now we are in the same state regardless of
        # whether the "-auto" flag was passed or not. The "Units" file must
        # be read in by "add_units_from_file" in method "populate" below.
        self.units = {}

    def populate(self):
        # main stuff
        self.check_empty()
        sim_unit = [gvars.simulation_dir]

        # read Units file and add units found there
        self.add_units_from_file()
        if not gvars.auto:
            # check if a Simulation is already present
            names = [x for x in self.get_list() if \
                     any([x.startswith(path) for path in gvars.simulations_path])]
            try:
                # remove unit if present
                names.remove(sim_unit)
            except ValueError:
                pass

            if names:
                gvars.out.put("\nERROR: Multiple Simulation Units Found", setup_globals.ERROR)
                gvars.out.push()
                gvars.out.put("You tried to setup {0} as well as ".format(sim_unit),
                              setup_globals.ERROR)
                gvars.out.put("use your previous (or specified) Units file.", setup_globals.ERROR)
                gvars.out.put("Your units file already has "
                              "{0}. Your choices are \n".format(names[0]), setup_globals.ERROR)
                gvars.out.put("* Use the -auto option and I will use my crystal ball",
                              setup_globals.ERROR)
                gvars.out.put("* Explicitly specify the units file using -unitsfile option",
                              setup_globals.ERROR)
                gvars.out.put("* Change the problem you want to setup", setup_globals.ERROR)
                gvars.out.pop()
                raise SetupError("")
            else:
                self.add_units(sim_unit)

        gvars.out.put('checking for default sub-units in included units', setup_globals.INFO)

        self.add_default_units()
        self.add_parent_units()
        self.add_kernel_units()

        gvars.out.put('checking for exclusivity and conflicts', setup_globals.INFO)
        self.check_exclusivity()
        self.check_conflicts()

        gvars.out.put('looking for paths of non-included units', setup_globals.INFO)
        self.check_top_units()
        self.add_top_units()

        self.add_units(sim_unit)

        gvars.out.put('checking requirements', setup_globals.INFO)
        self.check_requirements()

        self.check_suggest()

        # remove specified units and children (USER BEWARE no CHECK performed on these)
        self.remove_subtrees(gvars.killUnits.keys())

    def get_rp_info(self, max_plot_vars, **kw):
        rp_info = RPInfo(**kw)
        for (unitname, unit) in self.units.items():
            for (rpname, (rptype, rpvalue, rpconst, rprange)) in unit.parameter.items():
                # try rpname, then rpname_parameter
                rpcomment = unit.d.get(rpname, unit.d.get("{0}_parameter".format(rpname), ""))
                rp_info.addRP(rpname, type=rptype, value=rpvalue, const=rpconst, range=rprange,
                             location=unitname, comment=rpcomment)
            # now for all documentation without corresponding parameter declarations
            # which start with __
            for (key, value) in unit.d.items():
                if key[:2] == '__' and key not in unit.parameter:
                    # a type of DOC is for documentation only
                    rp_info.addRP(key, type='DOC', location=unitname, comment=value)

        # Generate enough plot_var_N names for all variables in the simulation.
        if type(max_plot_vars) == type(1):
            numRegularPlotVars = max_plot_vars
            for i in range(1, numRegularPlotVars + 1):
                rp_info.addRP('plot_var_{0}'.format(i), type="STRING", value='"none"',
                             location="IO/IOMain", comment="(automatically generated by setup)")
        else:
            rp_info.addRP('plot_var_{0}'.format(max_plot_vars), type="STRING", value='"none"',
                         location="IO/IOMain", comment="(automatically generated by setup)")
        return rp_info

    def get_var_info(self, **kw):
        var_info = VarInfo(**kw)
        for (unitname, unit) in self.units.items():
            # Add all variables
            for (var, (attr, eos1, eos2)) in unit.variable.items():
                var = var.lower()
                varcomment = unit.d.get("{0}_variable".format(var), "")
                var_info.addVar(var, type="Variable", location=unitname,
                                   attribs=attr, comment=varcomment)

            # Add all particle properties
            for (var, attr) in unit.particleprop.items():
                varcomment = unit.d.get("{0}_particleprop".format(var), "")
                var_info.addVar(var, type="Particle Property", location=unitname,
                                   attribs=[attr], comment=varcomment)

            # Add all particle types
            for var in unit.particletype.keys():
                varcomment = unit.d.get("{0}_particletype".format(var), "")
                var_info.addVar(var, type="Particle Type", location=unitname,
                                   attribs=[], comment=varcomment)

            # Add all particle property -> grid var maps
            for (var, attr) in unit.particlemap.items():
                var_info.addVar(var, type="Particle Property Map",
                                   location=unitname, attribs=[attr])

            # Add all species
            for var in unit.species.keys():
                varcomment = unit.d.get("{0}_species".format(var), "")
                var_info.addVar(var, type="Species", location=unitname,
                                   attribs=[], comment=varcomment)

            # Add all face variables
            for var in unit.facevar.keys():
                varcomment = unit.d.get("{0}_facevar".format(var), "")
                var_info.addVar(var, type="Face Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            # Add all scratch variables
            for var in unit.scratchvar.keys():
                varcomment = unit.d.get("{0}_scratchvar".format(var), "")
                var_info.addVar(var, type="Scratch Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            for var in unit.scratchcentervar.keys():
                varcomment = unit.d.get("{0}_scratchcentervar".format(var), "")
                var_info.addVar(var, type="Scratch Center Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            for var in unit.scratchfacexvar.keys():
                varcomment = unit.d.get("{0}_scratchfacexvar".format(var), "")
                var_info.addVar(var, type="Scratch Face-X Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            for var in unit.scratchfaceyvar.keys():
                varcomment = unit.d.get("{0}_scratchfaceyvar".format(var), "")
                var_info.addVar(var, type="Scratch Face-Y Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            for var in unit.scratchfacezvar.keys():
                varcomment = unit.d.get("{0}_scratchfacezvar".format(var), "")
                var_info.addVar(var, type="Scratch Face-Z Variable", location=unitname,
                                   attribs=[], comment=varcomment)

            # Add all Mass Scalars
            for var in unit.mass_scalar.keys():
                varcomment = unit.d.get("{0}_mass_scalar".format(var), "")
                var_info.addVar(var, type="Mass Scalar", location=unitname,
                                   attribs=[], comment=varcomment)
            # Add all Fluxes
            for var in unit.flux.keys():
                varcomment = unit.d.get("{0}_flux".format(var), "")
                var_info.addVar(var, type="Flux", location=unitname,
                                   attribs=[], comment=varcomment)
        return var_info

    def adjust_opts(self):
        """set default values for options possibly based of list of units we are using"""
        gvars.out.put("Computing default values for options not specified on command line",
                      setup_globals.IMPINFO)

        # Find all Grid related units we have now
        gp = setup_globals.GRID_PREFIX
        if not gp.endswith(os.sep):
            gp = gp + os.sep
        gridList = [x[len(gp):] for x in self.get_list() if x.startswith(gp) and x[len(gp):]]

        # now gridList also includes non Grid subunits of "Grid/GridMain"
        # eliminate choices which do not occur in GRID_CHOICES
        grids = {}
        for x in gridList:
            # mark grids we have
            fp = x.split(os.sep)[0]
            if fp in setup_globals.GRID_CHOICES:
                grids[fp] = 1

        if not grids:
            gvars.out.put("No Grid specified/found", setup_globals.WARN)
            grid = ""
        elif len(grids) == 1:
            grid = grids.keys()[0]
            gvars.out.put("Using Grid {0}".format(grid), setup_globals.INFO)
        else:
            raise SetupError("Found multiple grids {0}!".format(str(grids.keys())))

        # Adjust maxblocks
        if gvars.maxblocks == None:
            # using UG
            if grid == 'UG':
                gvars.maxblocks = 1
            elif grid == 'paramesh':
                # using paramesh
                if gvars.dimension == 3:
                    gvars.maxblocks = 200
                else:
                    gvars.maxblocks = 1000
            else:
                # some other grid package
                if gvars.dimension == 3:
                    gvars.maxblocks = 200
                else:
                    gvars.maxblocks = 1000

        # check if we are trying nonfbs and not UG
        if not gvars.setup_vars.get("fixedBlockSize") and (grid != 'UG' and grid != 'Chombo'):
            gvars.out.put("Non-Fixed Block size works only with UG", setup_globals.ERROR)
            raise SetupError("Change your grid or switch to fixed block size")


    def write_setup_units_file(self):
        file = open(os.path.join(gvars.project_setup_dir, setup_globals.SETUP_UNITS_FILENAME), 'w')
        file.write("\n".join(self.get_list()))
        file.write("\n")
        file.close()

    def generate_setup_flash_units(self):
        fname = os.path.join(gvars.project_setup_dir, setup_globals.SETUP_FLASH_UNITS_FILENAME)
        tname = os.path.join(gvars.binDir, setup_globals.SETUP_FLASH_UNITS_TEMPLATE)

        tpl = Template(tname)
        tpl["unit_names"] = self.get_list()
        if len(tpl["unit_names"]) > 39:
            p_units = self.parent_units()
            c_units = [x for x in self.units.keys() if x not in p_units]
            if len(c_units) <= 39:
                # if compression actually got us into the allowed range...
                c_units.sort()
                gvars.out.put("Making list of units for %s more concise by removing parent units:\n   reduced %d -> %d unit names.".format(fname, len(tpl["unit_names"]), len(c_units)), setup_globals.INFO)
                tpl["unit_names"] = c_units
            else:
                # exactly one Capitalized component after zero or more small components?
                weedable = re.compile("^([a-z][^/]*/)*([A-Z][^/]*)$")
                alwayskeep = re.compile("Main$")  # never weed out a final *Main component
                c_units = [x for x in c_units if not weedable.match(x) or alwayskeep.match(x)]
                if len(c_units) <= 39:
                    # if more compression now got us into the allowed range...
                    c_units.sort()
                    gvars.out.put("Making list of units for %s more concise by removing parent and stub-only units:\n   reduced %d -> %d unit names.".format(fname, len(tpl["unit_names"]), len(c_units)), setup_globals.INFO)
                    tpl["unit_names"] = c_units
                else:
                    # exactly one Capitalized component after zero or more
                    # small components floowed by "/localAPI"?
                    weedable = re.compile("^([a-z][^/]*/)*([A-Z][^/]*)/localAPI$")
                    # never weed out a final *Main component
                    alwayskeep = re.compile("Main/localAPI$")
                    c_units = [x for x in c_units if not weedable.match(x) or alwayskeep.match(x)]
                    if len(c_units) <= 39:
                        # if more compression now got us into the allowed range...
                        c_units.sort()
                        gvars.out.put("Making list of units for %s more concise by removing localAPI directories:\n   reduced %d -> %d unit names.".format(fname, len(tpl["unit_names"]), len(c_units)), setup_globals.INFO)
                        tpl["unit_names"] = c_units
        tpl.generate(fname)

    def create_makefiles(self):
        """ Makefiles: one for each toplevel Unit. Makefiles for subunits get
        appended to the top one."""

        tname = os.path.join(gvars.binDir, setup_globals.MAKEFILE_STUB_TEMPLATE)
        tpl = Template(tname)
        for unitname in self.get_list():
            lowest_base = get_lowest_base(unitname)

            # if we find a unit (ie leaf is a capital
            # otherwise it's just an organizational directory
            # the exception is "flashUtilities" which is not a
            # Unit.  Look for it specifically
            if not lowest_base[0].isupper():
                if lowest_base.find("flashUtilities") >= 0:
                    lowest_base = "flashUtilities"
                else:
                    continue

            source = os.path.join(self.units[unitname].unitpath, "Makefile")
            target = 'Makefile.' + lowest_base
            target = os.path.join(gvars.project_setup_dir, target)

            if os.path.isfile(source):
                file = open(target, 'a')
                file.write(open(source).read())
                file.close()
            elif unitname == lowest_base:
                tpl["unitname"] = unitname
                tpl.generate(target)


#
# Internally called code
#
def get_lowest_base(base):
    """Return the part of name closest to ROOT dir, which starts with capital letter."""
    parts = os.path.normpath(base).split(os.sep)
    #import pdb; pdb.set_trace()
    x = [p for p in parts if p[0].isupper()]
    try:
        return x[0]
    except IndexError:
        return base


