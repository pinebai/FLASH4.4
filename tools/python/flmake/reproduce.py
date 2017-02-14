# Python imports
import os
import sys
import shutil
import subprocess
import json
import uuid
import collections
import hashlib

# local imports
# These must be relative imports for flmake reproduce
from .. import FLASH_SRC_DIR, MPIRUN_CMD, _metadata
from ..utils import warning
from ..dsl import runtime_parameters
from . import logger
from . import setup_parse
from . import setup_globals
from .setup_globals import gvars
from .utils import desc_cmd_metadata, vc_checkout, vc_patch


def _checkout_and_patch(src, info, diff, cmd):
    targ = vc_checkout(src, info, cmd=cmd)
    if 0 < len(diff):
        vc_patch(targ, diff, info[0])
    return targ


def _run_prev_cmd(cmd, desc, srcdirs, prjdirs):
    # setup source dir
    src_diff =  desc[cmd]['source_diff'] if 'source_diff' in desc[cmd] else ""
    srcdir = _checkout_and_patch(FLASH_SRC_DIR, desc[cmd]['source_version'], src_diff, cmd)
    srcdirs.append(srcdir)
    md = dict(_metadata)
    md["FLASH_SRC_DIR"] = srcdir
    md["FLASH_CLEAN_SRC_DIR"] = os.path.join(srcdir, '.clean')
    md["MPIRUN_CMD"] = MPIRUN_CMD
    md["version"] = desc[cmd]['source_version']
    srctools = os.path.join(srcdir, 'tools')
    with open(os.path.join(srctools, 'python', 'metadata.json'), 'w') as f:
        json.dump(md, f, indent=2)
    sys.path.insert(0, srctools)

    # setup project dir
    prj_diff =  desc[cmd]['project_diff'] if 'project_diff' in desc[cmd] else ""
    prjdir = _checkout_and_patch(os.getcwd(), desc[cmd]['project_version'], prj_diff, cmd)
    prjdirs.append(prjdir)
    from python.flmake.setup_globals import gvars as prev_gvars
    _defs = prev_gvars._defaults
    prev_gvars.project_setup_dir = gvars.project_setup_dir
    prev_gvars.project_build_dir = gvars.project_build_dir
    prev_gvars.project_source_dir = os.path.join(prjdir, _defs['project_source_dir'])
    prev_gvars.project_simulations_dir = os.path.join(prjdir, _defs['project_simulations_dir'])
    rcfile = os.path.join(prjdir, 'flashrc.py')
    if 0 < len(desc[cmd]['rcfile']):
        with open(rcfile, 'w') as f:
            f.write(desc[cmd]['rcfile'])

    # run the command
    _mod_cmd(cmd, desc, srcdirs, prjdirs)
    from python.flmake.main import main as prev_main
    if cmd == "setup":
        # FIXME hack until setup_parse is refactored.
        sys.argv = desc[cmd]['command']
    prev_main(rcfile, desc[cmd]['command'][1:])

    # clear modules
    sys.path.remove(srctools)
    for mod in sys.modules.keys():
        if mod.startswith('python'):
            del sys.modules[mod]

    return srcdirs, prjdirs

def main(ns, rc):
    """Reproduces a FLASH run from a description."""
    gvars.init(FLASH_SRC_DIR)
    desc_filename = ns.desc

    # reporduce previous commands
    with open(desc_filename) as desc_file:
        desc = json.load(desc_file)

    srcdirs, prjdirs = [], []
    for cmd in ['setup', 'build', 'run']:
        if cmd in desc:
            srcdirs, prjdirs = _run_prev_cmd(cmd, desc, srcdirs, prjdirs)

    return 0


def usage():
    """Constructs usage info and exits."""
    msg = ("usage: flmake reproduce [options] <flash_descr>\n\n"
           "Reproduces a flash history by checking out previous\n"
           "versions of the code and re-executing them with the\n"
           "original values.  Note: it is highly advised to run\n"
           "'flmake clean 3' prior to reproduce.")
    return msg


#
# Modifications specific to commands
#


def _mod_cmd(cmd, desc, srcdirs, prjdirs):
    modfuncs = {'build': _mod_build, 'run': _mod_run}
    modfunc = modfuncs.get(cmd, lambda *a: None)
    modfunc(cmd, desc, srcdirs, prjdirs)


def _mod_build(cmd, desc, srcdirs, prjdirs):
    _relink_setup(gvars.project_setup_dir, srcdirs, prjdirs)


def _mod_run(cmd, desc, srcdirs, prjdirs):
    _relink_setup(gvars.project_setup_dir, srcdirs, prjdirs)


def _relink_setup(setupdir, srcdirs, prjdirs):
    allfiles = [os.path.join(setupdir, f) for f in os.listdir(setupdir)] 
    symlinks = [f for f in allfiles if os.path.islink(f)]
    srclinks = set([f for f in symlinks if os.readlink(f).startswith(srcdirs[-2])])
    prjlinks = set([f for f in symlinks if os.readlink(f).startswith(prjdirs[-2])])
    for f in symlinks:
        thedirs = srcdirs if f in srclinks else (prjdirs if f in prjlinks else None)
        oldtarg = os.readlink(f)
        newtarg = oldtarg.replace(*thedirs[-2:])
        os.remove(f)
        os.symlink(newtarg, f)
