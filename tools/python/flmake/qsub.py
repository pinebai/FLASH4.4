import os
import re
import sys
import shlex
import socket
import subprocess
import collections

if sys.version_info[0] == 3 or sys.version_info[:2] == (2, 7):
    import argparse
else:
    from .. import _argparse as argparse

# Relative imports needed!
from . import logger
from .setup_globals import gvars
from .utils import add_simple_opts
from .main import commands as main_commands, _make_argparser
from ..dsl import runtime_parameters

def _canonical_hostname(hostname):
    if hostname.endswith('eureka.alcf.anl.gov'):
        return 'eureka'
    elif hostname.endswith('beagle.ci.uchicago.edu'):
        return 'beagle'
    return hostname

HOSTNAME = _canonical_hostname(socket.getfqdn())
WHITESPACE = re.compile('\s')
QSUB_FILTER_OPTS = {
    'beagle': set(['t', 'n', 'np']),
    }

def qsub_command(arg_lists):
    """Get the FLASH qsub command qsub, based on a merger of argument lists."""
    # generate parser for qsub
    qsub_parser = argparse.ArgumentParser()
    for arg_list in arg_lists:
        add_simple_opts(qsub_parser, arg_list)

    # parse the options into the namespace, order matters!
    ns = argparse.Namespace()
    for arg_list in arg_lists:
        ns = qsub_parser.parse_args(arg_list, ns)

    # perform hostname specific modifications if needed
    execpath = qsub_hostname_mod(ns)
    if execpath is None:
        execpath = os.path.join(os.path.abspath(gvars.run_dir), 'flash4')

    # create the qsub command line string from the namespace
    opt_strings = dict([(act.dest, act.option_strings[0]) for act in \
                         qsub_parser._actions])
    opt_filter = QSUB_FILTER_OPTS.get(HOSTNAME, set())
    cmd = ['qsub']
    for opt, value in ns.__dict__.items():
        if opt in opt_filter:
            continue
        cmd.append(opt_strings[opt])
        if value is None:
            continue
        elif WHITESPACE.search(value):
            value = '"' + value + '"'
        cmd.append(value)
    cmd.append(execpath)
    return cmd, ns


def main(ns, rc):
    """submits job to scheduler with qsub."""
    # split command line arguments based on '--'
    dbl_dash_cnt = ns.options.count('--')
    if 0 == dbl_dash_cnt:
        run_args = ns.options[:]
        qsub_args = []
    elif 1 == dbl_dash_cnt:
        dbl_dash_ind = ns.options.index('--')
        run_args = ns.options[:dbl_dash_ind]
        qsub_args = ns.options[dbl_dash_ind+1:]

    # get the rc qsub
    qsub_rc = rc.get('qsub', [])
    if isinstance(qsub_rc, basestring):
        qsub_rc = shlex.split(qsub_rc)

    # run the run command in dry-run mode, using a dynamic import
    cmdmod, mainfunc = main_commands[ns.run_cmd]
    cmdmod = __import__(cmdmod, globals(), locals(), fromlist=[None])
    mainfunc = getattr(cmdmod, mainfunc)

    run_cmd_parser = _make_argparser()
    run_cmd_ns = run_cmd_parser.parse_args([ns.run_cmd] + run_args)
    run_cmd_ns.dry_run = True
    rtn = mainfunc(run_cmd_ns, rc)

    # log this attempted qsub run
    msg = ns.message
    if msg is None:
        msg = "scheduled run with: "
        msg += str(ns)
    logger.info(msg, "qsub", gvars.run_id, gvars.run_dir)

    # run qsub
    cmd, ns = qsub_command([qsub_rc, qsub_args])
    if hasattr(ns, 't'):
        qsub_walltime = float(ns.t)
        params = runtime_parameters.load(os.path.join(gvars.run_dir, 'flash.par'))
        params_wall = qsub_walltime - 15.0 if 45.0 <= qsub_walltime else qsub_walltime 
        params['wall_clock_time_limit'] = params_wall * 60
        runtime_parameters.dump(params, os.path.join(gvars.run_dir, 'flash.par'))
    rtn = subprocess.check_call(cmd, cwd=gvars.run_dir)
    return rtn


def usage():
    """Print usage info and exits."""
    msg = ("usage: flmake qsub <run-cmd> [-m MSG] [--] [<run-cmd options>]"
           " [--] [<qsub options>]\n\n"
           "Submits FLASH runs to the system scheduler using\n"
           "the 'qusub' utility. The user must specify which\n"
           "flmake run command (run, restart, etc) is desired.\n"
           "This run command is then executed in '--dry-run' mode\n"
           "with any other arguments & options that are needed.\n"
           "The qsub command may take optional arguments as well.\n"
           "These may be read from a string or list in flashrc.py and\n"
           "may optionally be given on the command line after '--',\n"
           "which signals the end of the run command options.\n\n")

    run_header = "QSUB HELP:\n"
    run_header += "-" * (len(run_header) - 1)
    msg += run_header + "\n"
    if hasattr(subprocess, 'check_output'):
        msg += subprocess.check_output(['qsub'], stderr=subprocess.STDOUT)[:-1]
    return msg


# Hostname modifications below
# This may get hairy!

def qsub_hostname_mod(ns):
    """Dispatch function for making hostname-specific modifications to how 
    qsub submits jobs.  This takes a namespace of all options and returns 
    the path to the executable to script we wish to run."""
    glbs = globals()
    hostfunc = HOSTNAME + "_mod"
    if hostfunc in glbs:
        f = glbs[hostfunc]
    elif hostfunc.replace('.', '_') in glbs:
        f = glbs[hostfunc.replace('.', '_')]
    else:
        f = lambda x: None
    execpath = f(ns) if callable(f) else None
    return execpath


def eureka_mod(ns):
    """Changes needed to make qsub run on Eureka."""
    runpath = os.path.abspath(gvars.run_dir)
    flshpath = os.path.join(runpath, 'flash4')
    execpath = os.path.join(runpath, 'run_script.sh')
    nodes = int(ns.n) if hasattr(ns, 'n') else 1
    nodes = int(ns.np) if hasattr(ns, 'np') else nodes
    procs = 8 * nodes  # Eureka has 8 procs per node 
    with open(os.path.join(gvars.project_setup_dir, 'Makefile.h')) as f:
        mpipath = [l for l in f if l.startswith("MPI_PATH")][0]
    mpipath = mpipath.partition('=')[-1].strip()
    withmx = 'mx' in os.path.split(mpipath)[-1]
    if withmx:
        script = ("#! /bin/sh\n"
                  "mpdboot -n {nodes} -f $COBALT_NODEFILE\n"
                  "mpiexec -n {procs} {flshpath}\n")
    else:
        script = ("#! /bin/sh\n"
                  "mpirun -machinefile $COBALT_NODEFILE -np {procs} {flshpath}\n")
    script = script.format(nodes=nodes, procs=procs, flshpath=flshpath)
    with open(execpath, 'w') as f:
        f.write(script)
    os.chmod(execpath, 0o755)
    return execpath


def beagle_mod(ns):
    """Changes needed to make qsub run on Beagle."""
    runpath = os.path.abspath(gvars.run_dir)
    flshpath = os.path.join('.', 'flash4')
    execpath = os.path.join(runpath, 'flash.pbs')
    nodes = int(ns.n) if hasattr(ns, 'n') else 1
    nodes = int(ns.np) if hasattr(ns, 'np') else nodes
    procs = 24 * nodes  # Beagle has 24 procs per node 
    walltime = float(ns.t)
    script = ("#!/bin/bash\n"
              "#PBS -N {run_dir}\n"
              "#PBS -l mppwidth={procs}\n"
              "#PBS -l mppnppn=24\n"
              "#PBS -l walltime={hours:02d}:{mins:02d}:{secs:02d}\n"
              "#PBS -V\n"
              "cd $PBS_O_WORKDIR\n"
              "flmake email --start\n"
              "aprun -B {flshpath} > stdout.txt\n"
              "flmake email --stop\n")
    script = script.format(nodes=nodes, procs=procs, flshpath=flshpath, 
                           run_dir=gvars.run_dir, hours=int(walltime)/60, 
                           mins=int(walltime)%60, secs=int((walltime%1)*60))
    with open(execpath, 'w') as f:
        f.write(script)
    os.chmod(execpath, 0o755)
    return execpath
