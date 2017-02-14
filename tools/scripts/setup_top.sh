#!/bin/bash
#
# Create a top-level directory in the current directory, suitable
# for running large Flash jobs with restarts by multiple people, and
# with provisions made for preserving run history.
#
# By Carlo Graziani, first draft 23 Jan 2007
#     Jan 2009 Klaus Weide  Minor shell syntax robustness improvement (for use on Intrepid)
#
# USAGE:
#
# setup_top.sh directory_name [options]
#
# options have the form "var=value".  They change script defaults.  Here
# are the defaults:
#

GID=$FLASH_GROUP_ID            # default set from the environment variable.
FLASH_PATH=$FLASH_SOURCE_TREE  # default set from the environment variable.
SITE=$FLASH_SITE_NAME          # default set from the environment variable.

#
[ $# -eq 0 -o "$1" = "--help" -o "$1" = "-h" ] && {

  echo "USAGE: setup_top.sh directory_name [options]"
  echo "See text of script for options"
  exit 1
}

DIRNAME=$1 ; shift

# Adjust defaults from command-line
while [ $# -gt 0 ] ; do eval $1 ; shift ; done

# Paranoia:
unalias -a
umask -S g+rwx

# Check that we actually have set $GID
[ "$GID" ] || {

  echo "No Group ID has been set.  Either set GID on the command line, or"
  echo "set the FLASH_GROUP_ID environment variable (why isn't it set in"
  echo "your shell setup files, apropos of nothing?)"
  exit 2
}

# ...and that we are members of this alleged group...
( groups | grep $GID > /dev/null ) || {

  echo "Regrettably, you do not appear to be a member of group $GID.  Much"
  echo "depends on your membership in said group, so either specify another"
  echo "group on the command-line, or acquire membership in $GID."
  exit 3
}

# ...and that the source tree is specified...
[ "$FLASH_PATH" ] || {

  echo "No path to the Flash source tree has been set.  Either set FLASH_PATH"
  echo "on the command line, or set the FLASH_SOURCE_TREE environment variable" 
  echo "(why isn't it set in your shell setup files, apropos of nothing?)"
  exit 4
}

# ... and that it exists...
[ -d "$FLASH_PATH" ] || {

  echo "The path $FLASH_PATH doesn't exist, or is not a directory, so it can't"
  echo "logically represent a Flash source code tree.  Please either check out"
  echo "a Flash tree at that path, or select another path."
  exit 5
}

# OK, let's set up.


mkdir $DIRNAME || {

  echo "Failed to create the directory $DIRNAME, a debugging opportunity."
  exit 6
}

chmod g+rws $DIRNAME || {

  echo "For some reason, it is not possible to set the directory permissions"
  echo "as required for consistent group access.  A debugging opportunity."
  exit 7
}

cd $DIRNAME

# Create the (empty) master operator notes file, as a hint.
touch MASTER_OPERATOR_NOTES

# Make a link to the source tree.
ln -s $FLASH_PATH flash_source

# See if we can locate a default batch submission script for this site.
if [ ! "$SITE" ]; then

  echo "No Site has been set.  Either set SITE on the command line, or"
  echo "set the FLASH_SITE_NAME environment variable (why isn't it set in"
  echo "your shell setup files, apropos of nothing?)"

elif [ ! -f flash_source/sites/$SITE/batch_submission_script ]; then

  echo "An appropriate batch submission script was not located in the source"
  echo "tree."

else

  cp flash_source/sites/$SITE/batch_submission_script .
  
fi

cat > OPERATOR_GUIDELINES << EOOG

************    OPERATOR BEST PRACTICES:


(1) Operators should ensure that they source the necessary
flash_environment script. These environment files should be in the flash
source tree, in the directory sites/<Current Site>.  They are called
"flash_environment.csh" and "flash_environment.sh".  Source whichever one
is appropriate to your shell.  These files will also append a path component
to \$PATH that will make available the scripts referenced below.

It is particularly important that every operator's file creation mask
should be set correctly -- "002" ("007" is OK too).  This will ensure that
all files created by operators are group-readable and group-writeable. 
Failure to do this will probably result in the creation of critical files
that nobody but the guilty operator may edit.  Type "umask" now, to see how
your umask is set.

Note that some sites (like Livermore, for example) set the umask to the
wrong value in _both_ the .cshrc and the .login file.  Livermore does this
using the local files .cshrc.blue, .login.blue (for up) and .cshrc.linux,
.login.linux (for alc), with the result that the umask is set correctly
in the .cshrc file, then immediately reset wrong by the .login file.  The
bottom line is, fix the umask in all your start-up scripts, if you find
that you keep creating files with the wrong permissions.

************ 

(2) Operators should create the top-level run directory using setup_top.sh.
The syntax is simply

  setup_top.sh <directory>

This will (a) create the directory, (b) set its group ownership to the
local Flash group ID, (c) ensure that it is group-readable, (d) ensure that
it is group-writeable, and (e) set it's SGID bit.  (e) will ensure
that all files created under it will also have the same group ownership,
and that all directories created under it will have the same attributes
set.  The script will then copy or link a few files in that directory.

The script (the source of this text, incidentally) has a few options that
can be controlled from the command line.  For example, to set the default
Flash source tree to something other than the project default, do

  setup_top.sh <directory> FLASH_PATH=/path/to/your/source/tree

************ 

(3) Operators should then copy a flash.par file into the newly-created
directory from somewhere (and possibly customize it).  That flash.par will
become the default flash.par that is copied into each run sub-directory.

************ 

(4) Operators should review the file "batch_submission_script" in this
directory, if one was found and imported, or else import an appropriate
one.  This is the script to be submitted to the local batch queue system.

************ 

(5) Operators should also edit the file MASTER_OPERATOR_NOTES.  Details
about the current series of computations -- what is to be accomplished, in
what priority, start-up parameters, flash.par entries worth scrutinizing
before each run, etc.  Anything that might serve as guidance to an
operator, particularly to one who isn't an expert, and can't get in touch
with one. Generally speaking, the contents of this file should be prepared
by team consensus well in advance of the run.

************ 

(6) The procedure for a normal start or re-start is to run setup_rundir.sh
from inside the directory created using setup_top.sh.  This will create a
new run directory.  The operator should then edit (or at least review) the
flash.par file in that directory as appropriate, and start the run.

(6a) The syntax for a new run -- NOT a restart -- is basically

  setup_rundir.sh

with no arguments.  This will create a new directory called (by default)
"rundir_XXX", where XXX are digits.  XXX is incremented each time the
script is run. The default flash.par (imported into the top-level directory
at step 3) is copied into the run directory after editing to ensure that 

  (a) restart = .false.
  (b) any plotFileNumber and ParticleFileNumber directives are commented
      out
  (c) the basenm supplied in the default flash.par WILL HAVE THE STRING
      "rundir_XXX_" PREPENDED TO IT to ensure a well-separated namespace
      for output files.

The default batch submission script is also copied into the new run
directory.

The default flash source tree is searched for a build directory called
"object", where the flash4 executable and any required data files are
expected.  Those files will be copied into the run directory.

The default flash source tree will also be queried (using svn) for version
info and modification status.  That information will also be placed in the
run directory.

The default assumptions of the script can be changed on the command line. 
For example, to use a different source tree from the default, and a
different name for the build directory,

  setup_rundir.sh FLASH_PATH=/path/to/your/source/tree \\
                  FLASH_BUILD_DIR=name_of_flash_build_dir

Consult the script text (or type "setup_rundir.sh --help") for other
options.

(6b) The procedure for a restart is to supply the script with BOTH the name
of a checkpoint file AND the name of a flash log file for the run that
produced that checkpoint file.  The syntax is:

  setup_rundir.sh RESTART_CHK=/path/to/checkpoint_file \\
                  RESTART_LOG=/path/to/flash_log_file

The behavior of the script is similar to the case of (6a), except that
in the flash.par that is copied to the run directory,

 (a) restart = .true.
 (b) the checkpointFileNumber is adjusted to suit the supplied checkpoint
     file
 (c) the plotFileNumber and ParticleFileNumber will be adjusted so
     as to be incremented by one with respect to the respective numbers
     immediately preceding the checkpoint file, as determined by parsing
     the supplied flash log file.

The supplied checkpoint file is symlinked with a link name chosen for
consistency with the value of "basenm".

************ 

(7) By default, it is intended that operators should do flash setups and
makes in the default Flash source directory.  The file 'flash_source' is
actually a symlink to that directory.  If this is done, then
setup_rundir.sh should correctly pick up the flash4 executable, the
required data files, the configuration history, and the subversion info
necessary for any required post-mortems.

If you make local modifications to that tree before you compile, that's OK,
since the changes will be noticed by subversion and recorded in the run
directory.

If you want to use your own tree for all runs, rather than the default
project tree, either tell setup_top.sh's command line where that tree is,
or modify the symlink "flash_source" to point to your tree.

To use a custom tree for a single run, tell setup_rundir.sh where to find
it on the command line using FLASH_PATH=/path/to/your/source/tree.

Please, please, do not just copy in a random flash4 executable from some
unrecorded tree with undocumented and uncommitted source modifications.  If
something goes wrong -- or worse, if we notice three weeks later that
something was wrong -- we probably wouldn't be able to reconstruct what
exactly was in that executable, which might make identifying the source of
the problem nearly impossible. Not keeping a record of the contents of the
executable is like performing lab work without keeping a lab notebook.
Please don't do that.

************ 

(8) Operators should avoid restarts in the same directory as a previous
non-trivial run.  A "non-trivial" run is one in which Flash actually
computed something worth keeping, or exited immediately with an interesting
error that will require later analysis.

This is not to say that one might never wish to restart in the current
directory.  If, for example, a parameter error were detected soon after the
run started, and the operator decides to kill the job and restart with the
correct parameters, it is not in principle obnoxious to recycle the current
run directory.  It would be a good idea to delete any output from the
aborted run before restarting, as we appear to have architecture-related
file clobbering issues due to at least one implementation of MPI or HDF (on
up) not respecting the file creation mask.  This issue can cause a run to
fail when it attempts to over-write an existing checkpoint or plot or
particle file.

************ 

(9) Operators should start Bob's Watcher Script when a run is on its way.
This may possibly be started from the batch submission script on some
systems.

************ 

(10) Each run directory will receive an (empty) file called
OPERATOR_NOTES.  The OPERATOR_NOTES file is for the sake of preserving
information about what the restart was about.  A few words should be
sufficient, e.g. "Previous run crashed due to excessive block counts, so we
re-run with nrefine_max throttled back".  The idea here is not a detailed
description, but rather enough information to help someone digging around
later to figure out what this run was for without having to inspect a lot
of Flash input and output.

EOOG
