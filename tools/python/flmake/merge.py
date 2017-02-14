import os
import shutil
import json

# relative imports needed!
from . import logger
from .run import init_rundir
from .setup_globals import gvars
from .utils import desc_cmd_metadata

USAGE = ("Merges a run directory with all previous\n"
         "runs in its history.  If a target directory\n"
         "is provided, will perform merge in this dir.\n"
         "Does nothing if a top-level run (no history).\n"
         "Merges get new ids and empty histories.\n\n"
         "usage: flmake merge <leaf-dir|id> [<target-dir>]")

def main(ns, rc):
    """Merges run with its history."""
    init_rundir()

    # get run dirs that 
    subdirs = [f for f in os.listdir('.') if os.path.isdir(f)]
    rundirs = [d for d in subdirs if os.path.exists(os.path.join(d, gvars.desc_filename))]

    # grab the hashes, diriectory names, and descriptions
    id_dir = {}
    id_desc = {}
    for rundir in rundirs:
        with open(os.path.join(rundir, gvars.desc_filename)) as desc_file:
            desc = json.load(desc_file)

        if 'run' not in desc:
            continue

        id_dir[desc['run']['id']] = rundir
        id_desc[desc['run']['id']] = desc

    # Get leaf id and dir
    leaf_id = ns.leaf if ns.leaf in id_dir else \
              [id for id, d in id_dir.items() if os.path.samefile(ns.leaf, d)][0]
    leaf_dir = id_dir[leaf_id]

    # Get target dir
    target_id = gvars.run_id
    target_dir = ns.target if ns.target is not None else gvars.run_dir
    if os.path.exists(target_dir):
        raise SystemExit("target directory {0} already exists.".format(target_dir))
    os.mkdir(target_dir)

    # Copy over leaf files to target dir
    try:
        for f in sorted(os.listdir(leaf_dir)):
            src = os.path.join(leaf_dir, f)
            dst = os.path.join(target_dir, f)
            print "{0} -> {1}".format(src, dst)
            shutil.copy2(src, dst)
    except shutil.Error as e:
        shutil.rmtree(target_dir)
        raise e

    # Get history ids from leaf
    with open(os.path.join(leaf_dir, gvars.desc_filename)) as leaf_desc_file:
        leaf_desc = json.load(leaf_desc_file)
    history = leaf_desc['run']['history'] + [leaf_desc['run']['id']]

    # Copy over metadata from all previous runs in history
    metafiles = set([gvars.desc_filename, 'flash.par'])
    rename = lambda i, f: '.'.join([f.rpartition('.')[0], i, f.rpartition('.')[2]])
    srcdst = [(os.path.join(id_dir[id], mf), os.path.join(target_dir, rename(id, mf))) \
              for id in history for mf in metafiles]
    for src, dst in srcdst:
        print "{0} -> {1}".format(src, dst)
        shutil.copy2(src, dst)

    # Open previous description
    with open(os.path.join(target_dir, gvars.desc_filename)) as desc_file:
        desc = json.load(desc_file)

    # write new flash description 
    prev_desc_run = desc['run']
    desc['run'] = desc_cmd_metadata()
    desc_run = desc['run']

    desc_run['id'] = gvars.run_id
    desc_run['history'] = []

    if 'merge' not in desc:
        desc['merge'] = {}

    if 'history' not in desc['merge']:
        desc['merge']['history'] = []

    desc['merge']['history'] += history

    with open(os.path.join(target_dir, gvars.desc_filename), 'w') as f:
        json.dump(desc, f, indent=2)

    # log this merge
    msg = ns.message
    if msg is None:
        msg = "merged {0} into {1} ({2})".format(", ".join(history), target_id, target_dir)
    logger.info(msg, "merge", target_id, target_dir)
