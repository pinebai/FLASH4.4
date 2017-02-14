# Python imports
import collections

# local imports
# relative imports needed
from .run import main as run_main


def main(ns, rc):
    """Sweeps over a suite of runs."""
    # ensure sweeping is list-like.
    rcsweep = rc['sweep']
    if callable(rcsweep):
        rcsweep = rcsweep()
    assert isinstance(rcsweep, collections.Sequence)

    # run the runs
    for s in rcsweep:
        rcs = dict(rc.items())
        if isinstance(s, collections.Mapping):
            ns.target = None
            rcs['parameters'] = s
        elif isinstance(s, collections.Sequence) and 2 == len(s):
            ns.target = s[0]
            rcs['parameters'] = s[1]
        rtn = run_main(ns, rcs)
    return rtn


USAGE = """usage: flmake sweep [--dry-run] [options]

Runs a suite of parameters given by the 'sweep' variable
in the run control file.  'sweep' may either be a list 
or a function which returns a list.  This list may contain
either parameter dictionaries (that the run command may
accept) or tuples of the run name and the parameter dict.
"""
