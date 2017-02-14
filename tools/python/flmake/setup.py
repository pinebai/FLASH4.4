# Python imports
import os
import sys
import re
import shutil
import types
import glob
import socket
import math
import json
import warnings

# local imports
# relative imports needed!
from .. import FLASH_SRC_DIR, _metadata
from ..utils import warning, message
from . import setup_parse
from . import setup_globals
from .setup_globals import gvars, SetupError
from .utils import dir_glob, determine_machine, search_paths, get_rel_path, \
     desc_cmd_metadata
from .lib_union import LibUnion
from .link_file_list import LinkFileList
from .setup_configuration import ConfigurationList
from .gen_files import generate_flash_defines, generate_makefile, \
     write_simulation_files, generate_pm3rp_defines, generate_buildstamp_generator

RUNTIME_FILES = ["flash.par", setup_globals.PM3_RP_FILENAME]

def base_setup_main(ns, rc):
    """Does the heavy lifting for the setup command."""
    global RUNTIME_FILES

    # this is not a threaded application
    sys.setcheckinterval(10000)
    cwdir = os.getcwd()

    # setup.py is FLASH_HOME/bin/setup.py
    # FLASH_HOME is the parent of directory having setup.py
    if FLASH_SRC_DIR:
        # New behaviour
        flash_home = FLASH_SRC_DIR
    else:
        # Old behavior
        flash_home = os.path.dirname(os.path.dirname(os.path.abspath(sys.argv[0])))

    # Initialize other related directory names
    gvars.init(flash_home)

    # parse the command line arguments. The effect of calling
    # this method is to set a lot of values in 'gvars'
    arg_lists = [rc['setup'], ns.options] if 'setup' in rc else [ns.options]
    setup_parse.parse(arg_lists)

    # Proceed with other stuff
    machdir = determine_machine()

    # A class which encapsulates operations of unit collections
    warnings.simplefilter('ignore')  # suppress the warnings in the config files
    unit_list = ConfigurationList()

    # Add local setup dir
    lsetdir = os.path.join(cwdir, gvars.project_setup_dir)
    if not os.path.isdir(lsetdir):
        os.mkdir(lsetdir)
    gvars.project_setup_dir = lsetdir

    # copy over the user specified units file
    if gvars.unitsfile:
        sname = os.path.join(gvars.flash_src_dir, gvars.unitsfile)
        tname = os.path.join(gvars.project_setup_dir, setup_globals.UNITS_FILENAME)
        shutil.copy(sname, tname)

    # generate and write an automatic "Units" file
    if gvars.auto:
        unit_list.generate_units_file()

    # read units file, satisfy requirement and all the good stuf
    unit_list.populate()

    # adjust values of command line options based on units used
    # this is the only place the setup script uses its knowledge
    # about Flash Units and what they do
    unit_list.adjust_opts()
    setup_parse.final()

    # create the object which does stuff relating to linking files
    link_list = LinkFileList(lsetdir)

    # Combine info from all Config files
    config_info = unit_list.get_config_info()

    # get Runtime Paramas info
    rp_info = unit_list.get_rp_info(config_info.max_plot_vars)

    # get info regarding all variables
    var_info = unit_list.get_var_info()

    # Add library options to config_info
    for key in gvars.with_libraries.keys():
        if key not in config_info.library:
            config_info.library[key] = []
        config_info.library[key].append(gvars.with_libraries[key])

    # Now handle libraries depending on other libraries
    # ConfigInfo.lib_config_info is a dictionary mapping
    # libraries to their config information
    # also CLASS.libOrder is a list of pairs (libname, args)
    # This order is important due to linker invocation order
    config_info.lib_config_info = LibUnion(config_info.library)

    # write out lib data to setup_libraries
    config_info.lib_config_info.write_libraries()

    # only now is noClobber flag stable
    link_list.clean_setup_dir()

    # remove the success file
    try:
        os.unlink(setup_globals.SUCCESS_FILENAME)
    except OSError:
        pass

    # write Fortran code for RP support
    rp_info.write_code(config_info)

    # find files which should not be linked
    link_list.get_dont_link_list(unit_list.get_list(), config_info.linkif)

    # get_link_order does the fancy sorting of unitnames
    for unitname in unit_list.get_link_order():
        link_list.link_files(unit_list.units[unitname].unitpath)
    link_list.link_files(os.path.join(gvars.flash_src_dir, 'sites'))

    # Link in the right version of the make file
    # attempt to find Makefile in FLASH root directory, otherwise link make file from sites
    makefile_path_root = os.path.join(flash_home, "Makefile.h" + gvars.makefileext)
    makefile_path_sites = os.path.join(machdir, "Makefile.h" + gvars.makefileext)

    if os.path.isfile(makefile_path_root):
        link_list.add_link(makefile_path_root, os.path.join(gvars.project_setup_dir, "Makefile.h"))
        print "Using Makefile.h: " + makefile_path_root
    else:
        link_list.add_link(makefile_path_sites, os.path.join(gvars.project_setup_dir, "Makefile.h"))
    print "Using Makefile.h: " + makefile_path_sites

    # functions in the simulations dir override earlier instances
    sim_dir = search_paths(gvars.simulations_path, gvars.simulation_name)
    sim_dir = os.path.abspath(sim_dir)
    sim_location = os.path.split(sim_dir)[0]
    link_list.link_files(sim_location)
    link_list.link_files(sim_dir)

    # functions from LINKIF statements override everything else,
    # assuming their conditions apply.
    link_list.do_linkif_overrides(unit_list.get_list(), config_info.linkif)

    # now is when we do the real link/copying
    link_list.really_link()

    #
    # Makefiles
    #
    test_src = os.path.join(gvars.flash_src_dir, 'tools', 'scripts', 'testing',
                            'precision_test', 'precision_test.F90')

    if os.path.isfile(test_src):
        shutil.copy(test_src, gvars.project_setup_dir)

    # FIXME merge with create_makefiles
    gvars.out.put('creating Makefiles for all units', setup_globals.INFO)
    unit_list.create_makefiles()

    gvars.out.put('generating buildstamp generator', setup_globals.INFO)
    generate_buildstamp_generator()

    gvars.out.put('copying release accessor function Makefile', setup_globals.INFO)
    shutil.copy(os.path.join(gvars.flash_src_dir, 'bin', 'make_release'), gvars.project_setup_dir)

    gvars.out.put('copying buildstats accessor function Makefile', setup_globals.INFO)
    shutil.copy(os.path.join(gvars.flash_src_dir, 'bin', 'make_bstats'), gvars.project_setup_dir)

    gvars.out.put('copying flashUnits accessor function Makefile', setup_globals.INFO)
    unit_list.generate_setup_flash_units()

    shutil.copy(os.path.join(gvars.flash_src_dir, 'bin', 'resetup'), gvars.project_setup_dir)

    gvars.out.put('copying Dependency generator', setup_globals.INFO)
    flist = ["setup_depends.py", "setup_addcdepends.py", "setup_reorder.py", "reorder.tpl"]
    for f in flist:
        shutil.copy(os.path.join(gvars.flash_src_dir, 'bin', f), gvars.project_setup_dir)

    gvars.out.put('generating Makefile', setup_globals.IMPINFO)

    # writes the Makefiles, Flash.h file, and simulation files.
    generate_makefile(config_info, machdir)
    generate_flash_defines(config_info)
    write_simulation_files(config_info)

    if (gvars.setup_vars.get("ParameshLibraryMode") or 
       ((not gvars.setup_vars.get("Grid") or (gvars.setup_vars.get("Grid").upper()=="PM4DEV")) 
        and ('USE_AMR_RUNTIME_PARAMETERS_FILE' in config_info.ppdefines or
        ('-DUSE_AMR_RUNTIME_PARAMETERS_FILE' in gvars.definesNames)) )):
        generate_pm3rp_defines(config_info)

    # Copy/link datafiles to object directory
    if gvars.portable:
        runtime_link_cmd = lambda real, link: shutil.copy2(real, link)
    else:
        cwd = os.getcwd()
        runtime_link_cmd = lambda real, link: os.symlink(os.path.abspath(real), link)

    # Fill runtuime files
    for wildcard in config_info.datafiles:
        if wildcard.startswith(gvars.project_simulations_dir):
            RUNTIME_FILES.extend([os.path.abspath(f) for f in glob.glob(wildcard)])
        else:
            RUNTIME_FILES.extend(glob.glob(os.path.join(gvars.source_dir, wildcard)))
    RUNTIME_FILES = [f for f in RUNTIME_FILES if not f.endswith('.tpl')]

    for f in RUNTIME_FILES:
        src_file = os.path.join(sim_dir, f)
        dst_file = os.path.join(gvars.project_setup_dir, os.path.basename(f))

        # special case certain entries
        if f == 'flash.par' and gvars.parfile:
            src_file = os.path.join(sim_dir, gvars.parfile)
            gvars.out.put("Copied {0} as flash.par".format(gvars.parfile), 
                          setup_globals.IMPINFO)
        elif f == setup_globals.PM3_RP_FILENAME:
            continue

        # link files
        if os.path.isfile(src_file) and not os.path.isfile(dst_file):
            runtime_link_cmd(src_file, dst_file)


    if RUNTIME_FILES:
        gvars.out.put('Copying data files: {0} copied'.format(len(RUNTIME_FILES)),
                      setup_globals.IMPINFO)

    setup_parse.write_cmd_line()
    unit_list.write_setup_units_file()
    var_info.write_var_info()

    # rp stuff
    rp_info.write_rp_info()
    rp_info.write_default_par()

    print message('SUCCESS')


def setup_desc(ns, rc):
    # Writes the flash description file for setup.
    desc = {'setup': desc_cmd_metadata()}
    desc_setup = desc['setup']
    desc_setup['datafiles'] = [os.path.basename(f) for f in RUNTIME_FILES]
    desc_setup['reproducible'] = len(desc_setup['source_version']) != 0 and \
                                 len(desc_setup['project_version']) != 0
    if not desc_setup['reproducible']:
        print warning("Irreproducible: FLASH source and/or project dirs not "
                      "under version control!")
    with open(gvars.desc_filename, 'w') as f:
        json.dump(desc, f, indent=2)


def main(ns, rc):
    """Initializes a new simulation."""
    base_setup_main(ns, rc)
    setup_desc(ns, rc)


def script_main(ns, rc):
    """A version of setup main that closely mimics the the old setup script."""
    base_setup_main(ns, rc)
    try:
        #base_setup_main(ns, rc)
        pass
    except SetupError, inst:
        if inst.args:
            print inst.args[0]
        sys.exit(1)
    except KeyboardInterrupt:
        print '\nuser abort'
        sys.exit(2)
    except Exception as e:
        if hasattr(e, 'dumbuser'):
            print '\nUser Error!!!\n{0}'.format(e.message)
        else:
            print '\nA setup internal error has occured, if possible please email the following'
            print 'debugging info to flash-bugs@flash.uchicago.edu'
            print 'Arguments:', sys.argv
            print 'Python Version: {0}.{1}.{2}'.format(*sys.version_info[:3])
            print 'Platform Details: {0}'.format(sys.platform)
            raise e


if __name__ == '__main__':
    main([], {}, "")
