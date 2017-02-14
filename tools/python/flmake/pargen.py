import os
import json
import shutil

from ..dsl.runtime_parameters import load, par_form
    
USAGE = ("Prints runtime parameters using a run control file (e.g. flashrc.py)\n\n"
         "usage: flmake pargen flashrc.py")

def main(ns, rc):
    """Generate a runtime parameters file from a run control file"""
    srcfile = ns.rc

    # open run-control file, if present
    rc = {}
    if os.path.isfile(srcfile):
        execfile(srcfile, rc, rc)
    else:
        raise ValueError("Source must be a file")
    params = rc["parameters"]
    if callable(params): params = params()

    pardict = par_form(params)

    # convert to string
    paritems = pardict.items()
    paritems.sort()
    paritems = ["{0} = {1}".format(k, v) for k, v in paritems]
    rawfile = "\n".join(paritems)

    print rawfile
