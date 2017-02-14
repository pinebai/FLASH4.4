# relative imports needed
from . import setup_parse
from . import setup_globals
from . import build
from . import run
from . import restart
from . import clean
from . import metadata
from . import diffpar
from . import lsruns
from . import merge
from . import mv
from . import rm
from . import log
from . import qsub
from . import reproduce
from . import flmail
from . import pargen
from . import sweep

messages = {
    'setup': setup_parse.usage,
    'build': build.usage,
    'run': run.usage,
    'restart': restart.usage,
    'merge': merge.USAGE,
    'clean': clean.USAGE,
    'metadata': metadata.USAGE, 
    'diffpar': diffpar.USAGE, 
    'ls-runs': lsruns.USAGE, 
    'mv': mv.USAGE,
    'rm': rm.USAGE,
    'log': log.USAGE,
    'qsub': qsub.usage,
    'reproduce': reproduce.usage,
    'email': flmail.USAGE,
    'pargen': pargen.USAGE,
    'sweep': sweep.USAGE,
    }


def get_help(opts):
    """Gets the help message."""
    if 0 == len(opts) or opts[0] == 'help':
        from .main import commands
        commands = dict([(cmd, getattr(__import__(mod, globals(), locals(), fromlist=[None]), 
                            mainfunc)) for (cmd, (mod, mainfunc)) in commands.items()])
        msg = ("flmake: the FLASH workflow utility\n\n"
               "usage: flmake [-m MESSAGE] cmd [command options...]\n\n"
               "The following commands are available:\n\n  ")
        msg += "\n  ".join(["{0:<10}{1}".format(k, commands[k].__doc__) for k in sorted(commands)])
        msg += ("\n\nA typical workflow is as follows:\n\n"
                ">>> flmake setup -auto Sedov\n"
                ">>> flmake build\n"
                ">>> flmake run -n 1\n"
                ">>> flmake clean 3")
    elif hasattr(messages[opts[0]], '__call__'):
        try:
            msg = messages[opts[0]]()
        except setup_globals.SetupError:
            msg = None
    else:
        msg = messages[opts[0]]
    return msg


def main(ns, rc):
    """Prints help or usage messages for the various commands."""
    msg = get_help(ns.options)
    if msg:
        print msg
