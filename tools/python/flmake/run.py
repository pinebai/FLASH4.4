# Python imports
import os
import shutil
import subprocess
import json
import uuid
import collections
import hashlib

# local imports
# relative imports needed
from .. import FLASH_SRC_DIR, MPIRUN_CMD
from ..utils import warning, message
from . import logger
from . import setup_parse
from . import setup_globals
from .setup_globals import gvars
from .utils import desc_cmd_metadata
from ..dsl import runtime_parameters

CP_OBJ_FILES = ['flash4']

def init_rundir():
    """Initializes run_dir and run_id in gvars."""
    gvars.init(FLASH_SRC_DIR)
    #setup_parse.parse() # FIXME Eats mpirun opts

    # use default run id
    if not gvars.run_id:
        gvars.run_id = uuid.uuid4().hex[:8]

    # use default rundir
    if not gvars.run_dir:
        gvars.run_dir = gvars.run_dir_prefix + gvars.run_id


def main(ns, rc):
    """Executes flash in run directory."""
    # init directory
    if ns.target is not None:
        gvars.run_dir = ns.target
    init_rundir()
    if not os.path.isdir(gvars.run_dir):
        os.mkdir(gvars.run_dir)

    # Copy the RC file to the run directory:
    shutil.copy2(ns.rc, os.path.join(gvars.run_dir,"flashrc.py"))

    # init runtime files with build files.
    cp_run_files = dict([(f, os.path.join(gvars.project_build_dir, f)) for f in CP_OBJ_FILES])

    # add datafiles (from setup) to runtime files needed
    with open(gvars.desc_filename) as desc_file:
        desc = json.load(desc_file)
    from_setup = dict([(f, os.path.join(gvars.project_setup_dir, f)) for f in desc['setup']['datafiles']])
    cp_run_files.update(from_setup)

    # copy over runtime files
    for f in cp_run_files:
        if os.path.isfile(cp_run_files[f]):
            shutil.copy2(cp_run_files[f], os.path.join(gvars.run_dir, f))

    # apply run control
    rcparams = {}
    if 'parameters' in rc:
        # ensure parameters is dict-like.
        rcparams = rc['parameters']
        if callable(rcparams):
            rcparams = rcparams()
        assert isinstance(rcparams, collections.Mapping)

        # update the parfile if params non-empty
        if 0 < len(rcparams):
            parfile = os.path.join(gvars.run_dir, 'flash.par')
            params = runtime_parameters.load(parfile) if os.path.isfile(parfile) else {}
            params.update(rcparams)
            runtime_parameters.dump(params, parfile)

    # Writes the flash description file for running.
    desc['run'] = desc_cmd_metadata()
    desc_run = desc['run']
    desc_run['id'] = gvars.run_id
    desc_run['history'] = []
    flash_executable = os.path.join(gvars.run_dir, "flash4")
    with open(flash_executable, 'rb') as f:
        desc_run['flash_executable_hash'] = hashlib.sha1(f.read()).hexdigest()
    desc_run['flash_executable_mtime'] = os.path.getmtime(flash_executable)
    desc_run['reproducible'] = desc['build'] and \
        desc_run['flash_executable_hash'] == desc['build']['flash_executable_hash']
    if not desc_run['reproducible']:
        print message("Irreproducible: flash executable modifided between build and run!")
    with open(os.path.join(gvars.run_dir, gvars.desc_filename), 'w') as f:
        json.dump(desc, f, indent=2)

    # quit sanely, if desired
    if ns.dry_run:
        return 0

    # log this attempted run
    msg = ns.message
    if msg is None:
        msg = "running flash"
        msg += " with: " + str(rcparams) if 0 < len(rcparams) else ""
    logger.info(ns.message, "run", desc_run['id'], gvars.run_dir)

    # run flash
    rtn = 0
    nprocs = ['-n', ns.nprocs] if ns.nprocs is not None else []
    cmd = [MPIRUN_CMD] + nprocs + ns.options + ["./flash4"]
    try:
        rtn = subprocess.check_call(cmd, cwd=gvars.run_dir)
    except subprocess.CalledProcessError:
        pass
    finally:
        map(os.remove, [f for f in os.listdir('.') if 0 == os.path.getsize(f)])
    return rtn


def usage():
    """Print usage info and exits."""
    msg = ("usage: flmake run [--dry-run] [options]\n\n"
           "Runs the flash executable by calling out to\n"
           "the '{0}' utility. All options given to this\n"
           "command are transparently passed down.\n"
           "FLASH is executed in a new run directory which\n"
           "is given a unique id number upon each call.\n"
           "If '--dry-run' is present, no execution is done.\n\n").format(MPIRUN_CMD)

    run_header = "{0} HELP:\n".format(MPIRUN_CMD.upper())
    run_header += "-" * (len(run_header) - 1)
    msg += run_header + "\n"
    if hasattr(subprocess, 'check_output'):
        msg += subprocess.check_output([MPIRUN_CMD, '--help'], stderr=subprocess.STDOUT)[:-1]
    return msg
