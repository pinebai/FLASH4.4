import os
import json

# Relative imports needed
from .. import FLASH_SRC_DIR
from .setup_globals import gvars
from .utils import hash_list_to_dict, hash_dict_to_str

USAGE = ("Searches through local sub-directories\n"
         "and lists those which contain flash runs\n"
         "(and the run hashes) in a tree.\n\n"
         "usage: flmake ls-runs")

def main(ns, rc):
    """Lists flash run dirs as a tree."""
    gvars.init(FLASH_SRC_DIR)

    # get run dirs that 
    subdirs = [f for f in os.listdir('.') if os.path.isdir(f)]
    rundirs = [d for d in subdirs if os.path.exists(os.path.join(d, gvars.desc_filename))]

    # grab the appropriate data out of the description files
    labels = {}
    hash_lists = []
    for rundir in rundirs:
        with open(os.path.join(rundir, gvars.desc_filename)) as desc_file:
            desc = json.load(desc_file)

        if 'run' not in desc:
            continue

        labels[desc['run']['id']] = rundir + os.path.sep
        hl = desc['run']['history'] + [desc['run']['id']]
        hash_lists.append(hl)

    # Create graph stucture
    hash_dict = {}
    hash_lists.sort()
    for hash_list in hash_lists:
        hash_dict = hash_list_to_dict(hash_list, hash_dict)

    # Make into string, print, and return
    s = hash_dict_to_str(hash_dict, labels)[:-1]
    if 1 < len(s):
        print s
