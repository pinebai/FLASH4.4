import os
import sys
import subprocess

# relative imports needed
from ..dsl import runtime_parameters

USAGE = ("Converts two runtime parameter files\n"
         "to a cannonical form and displays the\n"
         "differnce.\n\n"
         "usage: flmake diffpar <parfile1> <parfile2>")

def _load_parfile(parfile):
    if not os.path.isfile(parfile):
        sys.exit("{0} is not a regular file".format(parfile))

    try:
        params = runtime_parameters.load(parfile)
    except ValueError as e:
        print e
        sys.exit("{0} has error(s); may not be a runtime parameter file".format(parfile))

    return params

def main(ns, rc):
    """Diffs two runtime parameter files."""
    # load files
    params1 = _load_parfile(ns.par1)
    params2 = _load_parfile(ns.par2)

    # define key sets
    keys1 = set(params1.keys())
    keys2 = set(params2.keys())

    inboth = keys1 & keys2
    onlyin1 = keys1 - keys2
    onlyin2 = keys2 - keys1

    # init output
    diffs = []

    # diff inboth
    for key in sorted(inboth):
        sameval = (params1[key] == params2[key]) 
        if sameval:
            continue

        if 0 == len(diffs):
            diffs.append("Parameters with different values:")
            diffs.append('-' * len(diffs[-1]))

        diffs.append("{0}: {1} != {2}".format(key, params1[key], params2[key]))

    # diff in 1
    if 0 < len(onlyin1): 
        if 0 < len(diffs):
            diffs.append('')
        diffs.append('Parameters only in {0}:'.format(ns.par1))
        diffs.append('-' * len(diffs[-1]))

    diffs += ["{0}: {1}".format(key, params1[key]) for key in sorted(onlyin1)]

    # diff in 2
    if 0 < len(onlyin2): 
        if 0 < len(diffs):
            diffs.append('')
        diffs.append('Parameters only in {0}:'.format(ns.par2))
        diffs.append('-' * len(diffs[-1]))

    diffs += ["{0}: {1}".format(key, params2[key]) for key in sorted(onlyin2)]

    # output string 
    if 0 < len(diffs):
        diffstr = "\n".join(diffs)
        print diffstr
        sys.exit(1)  # adhere to exit convention for diff
