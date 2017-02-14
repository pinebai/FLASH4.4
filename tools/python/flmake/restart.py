# Python imports
import os
import shutil
import subprocess
import json
import hashlib

# local imports
# relative imports needed
from .. import FLASH_SRC_DIR, MPIRUN_CMD
from ..utils import warning
from ..dsl import runtime_parameters
from . import logger
from . import setup_parse
from . import setup_globals
from .setup_globals import gvars
from .utils import desc_cmd_metadata
from .run import init_rundir


def main(ns, rc):
    """Executes flash in restart directory."""
    # init directory
    if ns.target is not None:
        gvars.run_dir = ns.target
    init_rundir()
    if not os.path.isdir(gvars.run_dir):
        os.mkdir(gvars.run_dir)

    # Set run id
    prev_run_dir = ns.prev_run

    # link files from previous run
    dont_relink = set([gvars.desc_filename, 'flash.par'])
    link_files = set(os.listdir(prev_run_dir)) - dont_relink
    for f in link_files:
        src = os.path.abspath(os.path.join(prev_run_dir, f))
        dst = os.path.abspath(os.path.join(gvars.run_dir, f))
        os.symlink(src, dst)

    # write new flash description 
    with open(os.path.join(prev_run_dir, gvars.desc_filename)) as desc_file:
        desc = json.load(desc_file)
    prev_desc_run = desc['run']
    desc['run'] = desc_cmd_metadata()
    desc_run = desc['run']
    desc_run['id'] = gvars.run_id
    desc_run['history'] = prev_desc_run['history'] + [prev_desc_run['id']]
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

    # write new parameter file based on previous one
    params = runtime_parameters.load(os.path.join(prev_run_dir, 'flash.par'))

    checkpoint_nums = [int(f.partition('_chk_')[2]) for f in link_files if '_chk_' in f]
    checkpoint_nums.sort()
    params['checkpointFileNumber'] = checkpoint_nums[-1] if 1 <= len(checkpoint_nums) else 0

    params['restart'] = (1 <= len(checkpoint_nums))  # only restart if there is a checkpoint

    plot_nums = [int(f.partition('_plt_cnt_')[2]) for f in link_files \
                             if ('_plt_cnt_' in f) and ('_forced_' not in f)]
    plot_nums.sort()
    params['plotFileNumber'] = plot_nums[-1] if 1 <= len(plot_nums) else 0

    forced_plot_nums = [int(f.partition('_plt_cnt_')[2]) for f in link_files \
                                    if ('_plt_cnt_' in f) and ('_forced_' in f)]
    forced_plot_nums.sort()
    params['forcedPlotfileNumber'] = forced_plot_nums[-1] if 1 <= len(forced_plot_nums) else 0

    # FIXME: need to add particle file restarts here

    runtime_parameters.dump(params, os.path.join(gvars.run_dir, 'flash.par'))

    # quit sanely, if desired
    if ns.dry_run:
        return 0

    # log this attempted restart
    msg = ns.message
    if msg is None:
        msg = "restarted {0}".format(prev_desc_run['id'])
    logger.info(ns.message, "restart", desc_run['id'], gvars.run_dir)

    # run flash
    nprocs = ['-n', ns.nprocs] if ns.nprocs is not None else []
    cmd = [MPIRUN_CMD] + nprocs + ns.options + ["./flash4"]
    try:
        rtn = subprocess.check_call(cmd, cwd=gvars.run_dir)
    finally:
        map(os.remove, [f for f in os.listdir('.') if 0 == os.path.getsize(f)])
    return rtn


def usage():
    """Print usage info and exits."""
    msg = ("usage: flmake restart <prevdir> [--dry-run] [options]\n\n"
           "Restarts the flash executable from a previous run.\n"
           "Flash is run with the '{0}' utility. All options \n"
           "given to this command are transparently passed down.\n"
           "Flash is executed in a new run directory which\n"
           "is given a unique id number upon each call.\n"
           "If '--dry-run' is present, no execution is done.\n\n").format(MPIRUN_CMD)

    run_header = "{0} HELP:\n".format(MPIRUN_CMD.upper())
    run_header += "-" * (len(run_header) - 1)
    msg += run_header + "\n"
    if hasattr(subprocess, 'check_output'):
        msg += subprocess.check_output([MPIRUN_CMD, '--help'], stderr=subprocess.STDOUT)[:-1]
    return msg

