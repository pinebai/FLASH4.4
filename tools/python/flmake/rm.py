import os
import json
import shutil

# relative imports needed
from .. import FLASH_SRC_DIR
from . import logger
from .setup_globals import gvars

USAGE = ("Deletes a flash run directory permanently.\n\n"
         "usage: flmake rm <rmdir>")

def main(ns, rc):
    """Deletes a local run dir."""
    gvars.init(FLASH_SRC_DIR)
    rmdir = ns.src

    # grab id, if possible
    id = None
    desc_filename = os.path.join(rmdir, gvars.desc_filename)
    if os.path.exists(desc_filename):
        with open(desc_filename) as desc_file:
            desc = json.load(desc_file)

        if 'run' in desc:
            id = desc['run']['id']

    # move the dir
    shutil.rmtree(rmdir)

    # Log the move
    msg = ns.message
    if msg is None:
        msg = "deleted {0}".format(rmdir)
    logger.info(msg, "rm", id, rmdir)
