#!/usr/bin/env python2.7

import os
import sys
import json
import subprocess 
from distutils.core import setup

USER = os.environ['USER']

if sys.version_info[0] == 2:
    if sys.version_info[1] >= 7:
        check_output = subprocess.check_output
    else:
        check_output = lambda *args: subprocess.Popen(*args, 
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()[0]

def make_metadata(path):
    """Build a metadata file as JSON."""
    srcdir = os.path.abspath(os.path.join(os.path.split(__file__)[0], '..'))
    md = {"FLASH_SRC_DIR": srcdir,
          "FLASH_CLEAN_SRC_DIR": os.path.join(srcdir, '.clean'),
          "MPIRUN_CMD": 'mpirun', 
          }

    # build version string
    ls = os.listdir('.')
    if '.git' in ls:
        rev = check_output(['git', 'rev-list', '-1', 'HEAD']).strip()
    elif '.svn' in ls:
        # FIXME to include svn rev number
        rev = ''
    else:
        rev = ''

    # munge version name
    release_path = os.path.join(srcdir, 'RELEASE')
    if os.path.exists(release_path):
        with open(release_path, 'r') as f:
            rel = f.read().strip()
    else:
        print "RELEASE file not found at " + release_path
        rel = os.path.split(srcdir)[-1]
    md['version'] = rel.replace('${SVNBRIEF}', rev).strip()

    # write the metadata file
    with open(path, 'w') as f:
        json.dump(md, f, indent=2)


bc_template = """\
# flamke bash completion

_flmake()
{{
    local cur prev cmds setup_args metadata_args cur_cmd cur_cmd_index
    COMPREPLY=()
    cur="${{COMP_WORDS[COMP_CWORD]}}"
    prev="${{COMP_WORDS[COMP_CWORD-1]}}"
    cmds="{cmds}"
    metadata_args="{metadata_args}"
    setup_args="{setup_args}"

    # find the location of the command depending on if there is a message
    if [[ "${{COMP_WORDS[1]}}" == "-m" ]] ; then
        cur_cmd_index=3
    else
        cur_cmd_index=1
    fi
    cur_cmd="${{COMP_WORDS[$cur_cmd_index]}}"

    # complete the command 
    if [[ $cur_cmd_index == 1 && $COMP_CWORD == 1 ]] ; then
        # match first command
        COMPREPLY=( $(compgen -W "${{cmds}} -m" -- ${{cur}}) )
    elif [[ $cur_cmd_index == 3 && $COMP_CWORD == 3 ]] ; then
        # match third command
        COMPREPLY=( $(compgen -W "${{cmds}}" -- ${{cur}}) )
    elif [[ $cur_cmd == "help" ]] ; then
        # match help argument to other commands
        COMPREPLY=( $(compgen -W "${{cmds}}" -- ${{cur}}) )
    elif [[ $cur_cmd == "metadata" &&  ${{cur}} == -* ]] ; then
        # match metadata
        COMPREPLY=( $(compgen -W "${{metadata_args}}" -- ${{cur}}) )
    elif [[ $cur_cmd == "setup" ]] ; then
        # match setup to args
        COMPREPLY=( $(compgen -W "${{setup_args}}" -- ${{cur}}) )
    elif [[ $cur_cmd == "diffpar" ]] ; then
        # match files for diff'ing parameters
        COMPREPLY=( $(compgen -f ${{cur}}) )
    elif [[ $cur_cmd == "restart" ]] ; then
        # match dirs for restarting
        COMPREPLY=( $(compgen -d ${{cur}}) )
    elif [[ $cur_cmd == "mv" ]] ; then
        # match dirs for moving
        COMPREPLY=( $(compgen -d ${{cur}}) )
    elif [[ $cur_cmd == "rm" ]] ; then
        # match dirs for deleting
        COMPREPLY=( $(compgen -d ${{cur}}) )
    elif [[ $cur_cmd == "log" && $prev == "log" ]] ; then
        # match log options
        COMPREPLY=( $(compgen -W "-n" -- ${{cur}}) )
    elif [[ $cur_cmd == "reproduce" ]] ; then
        # match files for diff'ing parameters
        COMPREPLY=( $(compgen -f ${{cur}}) )
    fi

    return 0
}}
complete -o filenames -F _flmake flmake
"""

def make_bash_completion():
    """Writes a bash completion file for flmake."""
    if USER == 'root':
        bcfile = '/etc/bash_completion.d/flmake'
    else:
        import flash.flmake.__init__
        bcfile = os.path.split(flash.flmake.__init__.__file__)[0]
        bcfile = os.path.join(bcfile, 'flmake')
        msg = ("\nFor bash-completion please add the following lines to your ~/.bashrc:\n"
               "# Enable completion for flmake\n"
               "if [ -f {bcfile} ] ; then \n"
               "    source {bcfile}\n"
               "fi")
        print msg.format(bcfile=bcfile)

    # Get commands automatically
    # This requiers that flash be installed
    from flash import FLASH_SRC_DIR
    from flash.flmake.main import commands
    from flash.flmake import setup_parse
    from flash.flmake.setup_globals import SHORTCUT_CHAR, gvars
    gvars.init(FLASH_SRC_DIR)
    del gvars.out

    # build the setup args
    setup_args = ['-' + a for a in setup_parse.WITH_ARGS]
    #setup_args += ['--' + a for a in setup_parse.WITH_ARGS]
    setup_args += ['-' + a for a in setup_parse.WITHOUT_ARGS]
    #setup_args += ['--' + a for a in setup_parse.WITHOUT_ARGS]
    setup_args += setup_parse.ADDL_DEF_ARGS
    setup_args += [SHORTCUT_CHAR + a for a in setup_parse.getShortcuts().keys()]
    for d in gvars.simulations_path:
        if not os.path.isdir(d):
            continue

        for sim in os.listdir(d):
            if not os.path.isdir(os.path.join(d, sim)) or sim.lower() == 'unittest':
                continue

            if sim[0].isupper():
                setup_args.append(sim)
            elif sim[0].islower():
                subsims = [os.path.join(sim, subsim) for subsim in os.listdir(os.path.join(d, sim))
                           if os.path.isdir(os.path.join(d, sim, subsim)) and subsim[0].isupper()]
                setup_args.extend(subsims)
                
    setup_args.sort()

    # write template
    bc_values = {'cmds': " ".join(commands.keys()), 
                 'metadata_args': '-e --edit',
                 'setup_args': " ".join(setup_args), 
                 }
    bc = bc_template.format(**bc_values)
    with open(bcfile, 'w') as f:
        f.write(bc)
    

def run_setup():
    """Runs the setup command"""
    # Create metadata file
    mdpath = "python/metadata.json"
    make_metadata(mdpath)

    # Actually run setup
    setup(name="flash",
        version='0.1',
        description='A modular, parallel multiphysics simulation code capable of handling general compressible flow problems found in many astrophysical environments.',
        author='The Flash Center',
        author_email='eder@flash.uchicago.edu',
        url='http://flash.uchicago.edu/',
        packages=['flash', 'flash.flmake', 'flash.dsl', 'flash.flashfile'],
        package_dir={'flash': 'python'},
        package_data={'flash': ['metadata.json']}, 
        scripts=['scripts/flmake', 
                 'scripts/consdat', 
                 'scripts/convertspect3d',
                 'scripts/flashtimes',
                 'scripts/symmetry',],
        )

    # clean up metadata file
    if os.path.exists(mdpath):
        os.remove(mdpath)

    # Install bash completion, must be done after setup call
    make_bash_completion()

    # remove empty flash.log
    if os.path.exists('flash.log') and 0 == os.path.getsize('flash.log'):
        os.remove('flash.log')


if __name__ == "__main__":
    run_setup()
