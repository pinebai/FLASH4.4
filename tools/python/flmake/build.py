# Python imports
import os
import subprocess
import shutil
import json
import hashlib

# local imports
from .. import FLASH_SRC_DIR
from ..utils import warning
from . import setup_parse
from .setup_globals import gvars
from .utils import desc_cmd_metadata, vc_diff


def build_main(ns, rc):
    """Does the heavy lifting for the build command."""
    # Some setup
    gvars.init(FLASH_SRC_DIR)
    initial_files = set(os.listdir(gvars.project_setup_dir))
    cmd = ['make', "-j", ns.j]

    # Build
    try:
        rtn = subprocess.check_call(cmd, cwd=gvars.project_setup_dir)
    finally:
        # Get source files after build and take diff
        after_files = set(os.listdir(gvars.project_setup_dir))
        new_files = after_files - initial_files

        # Move the new files to the object dir
        if not os.path.isdir(gvars.project_build_dir):
            os.mkdir(gvars.project_build_dir)

        for nf in new_files:
            src_file = os.path.join(gvars.project_setup_dir, nf)
            dst_file = os.path.join(gvars.project_build_dir, nf)
            os.rename(src_file, dst_file)

    # Kill flmake if make didn't build properly (after copying)
    assert rtn == 0

    # gets the flash description dictionary and updates it
    with open(gvars.desc_filename) as desc_file:
        desc = json.load(desc_file)
    desc['build'] = desc_cmd_metadata()
    desc_build = desc['build']
    flash_executable = os.path.join(gvars.project_build_dir, "flash4")
    with open(flash_executable, 'rb') as f:
        desc_build['flash_executable_hash'] = hashlib.sha1(f.read()).hexdigest()
    desc_build['flash_executable_mtime'] = os.path.getmtime(flash_executable)
    desc_build['reproducible'] = desc['setup']['reproducible'] and \
                                 len(desc_build['source_version']) != 0 and \
                                 len(desc_build['project_version']) != 0
    if not desc_build['reproducible']:
        print warning("Irreproducible: FLASH source and/or project dirs not "
                      "under version control!")
    with open(gvars.desc_filename, 'w') as f:
        json.dump(desc, f, indent=2)



def main(ns, rc):
    """Builds the flash executable."""
    build_main(ns, rc)


def usage():
    """Print usage info and exit."""
    msg = ("usage: flmake build [options]\n\n"
           "Builds the flash executable by calling out to\n"
           "the 'make' utility. All options given to this\n"
           "command are transparently passed down to make.\n"
           "After make executes, the results of flash \n"
           "compilation are moved to the object directory.\n\n")
    msg += "MAKE HELP:\n"
    msg += "----------\n"

    if hasattr(subprocess, 'echeck_output'):
        msg += subprocess.check_output(['make', '--help'])[:-1]
    return msg
