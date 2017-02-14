
import re, sys

def libinfo(relLibDir="",absLibDir="",buildFlag="",args="",macros=[]):
    ans = {}
    args = args.lower()

    if not args: # we want regular mpi
       return {"EXTERNAL":"mpi"}

    elif re.match("mpitrace", args):
       return {"EXTERNAL":"mpitrace"}

    elif re.match("mpihpm_smp", args):
       return {"EXTERNAL":"mpihpm_smp"}

    elif re.match("mpihpm", args):
       return {"EXTERNAL":"mpihpm"}

    elif re.match("dummy", args):
       return {"INTERNAL":"dmpi"}

    else: #unknown option 
       print >>sys.stderr, 'Unknown MPI Variant "%s"' % args
       return {}
