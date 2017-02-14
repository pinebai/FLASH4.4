import os
import json
import shutil

# Relative imports needed!
from .. import FLASH_SRC_DIR
from . import logger
from .setup_globals import gvars

USAGE = ("Moves a flash run local sub-directory\n"
         "from src to dst.  Useful for managing\n"
         "many runs.\n\n"
         "usage: flmake mv <src> <dst>")

def main(ns, rc):
    """Moves src run dir to dst dir."""
    gvars.init(FLASH_SRC_DIR)

    # grab id, if possible
    id = None
    desc_filename = os.path.join(ns.src, gvars.desc_filename)
    if os.path.exists(desc_filename):
        with open(desc_filename) as desc_file:
            desc = json.load(desc_file)

        if 'run' in desc:
            id = desc['run']['id']

    # move the dir
    shutil.move(ns.src, ns.dst)

    # Log the move
    msg = ns.message
    if msg is None:
        msg = "moved {0} -> {1}".format(ns.src, ns.dst)
    logger.info(msg, "mv", id, ns.dst)
