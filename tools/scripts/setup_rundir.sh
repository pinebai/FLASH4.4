#!/bin/bash
#
# Create a run directory with all necessary files to run Flash and
# preserve reconstructable history about the run.  This script should
# be run from inside a directory created by setup_top.sh.
#
# By Carlo Graziani, first draft 24 Jan 2007
#
# USAGE:
#
# setup_rundir [RESTART_CHK=/path/to/checkpoint_file RESTART_LOG=/path/to/flash_log_file] \
#              [FLASH_PATH=/path/to/your/source/tree] \
#              [FLASH_PAR=/path/to/flash.par] \
#              [FLASH_BUILD_DIR=name_of_flash_build_dir] \
#              [BATCH_SCRIPT=/path/to/submission/script] \
#              [RUNDIR_STEM=some_filename_stem] \
#              [EXISTING_RUNDIR=number]
#
# Defaults:

RESTART_CHK=                    # By default, not a restart.
RESTART_LOG=

FLASH_PATH=./flash_source       # Default set up by setup_topdir.sh

FLASH_PAR=./flash.par           # Operator should have provided this default

BATCH_SCRIPT=./batch_submission_script

FLASH_BUILD_DIR=object          # Where we'll look for flash4 and its files

RUNDIR_STEM=rundir_             # Run directories named "rundir_000", 
                                #    "rundir_001"...
#

CHK_SUFF=hdf5_chk_
PLT_SUFF=hdf5_plt_cnt_
PART_SUFF=hdf5_part_


# Paranoia:
unalias -a
umask -S g+rwx

[  $# -eq 1 ] && [ $1 = "--help" -o $1 = "-h" ] && {

  echo "USAGE: setup_rundir.sh [RESTART_CHK=/path/to/checkpoint_file \ "
  echo "                        RESTART_LOG=/path/to/flash_log_file] \ "
  echo "                       [FLASH_PATH=/path/to/your/source/tree] \ "
  echo "                       [FLASH_PAR=/path/to/flash.par] \ "
  echo "                       [BATCH_SCRIPT=/path/to/submission/script] \ "
  echo "                       [FLASH_BUILD_DIR=name_of_flash_build_dir] \ "
  echo "                       [RUNDIR_STEM=some_filename_stem]"
  echo "                       [EXISTING_RUNDIR=number]"
  echo ""
  echo "For defaults, consult the text of this script."
  exit 0
}

# Adjust defaults from command-line
while [ $# -gt 0 ] ; do eval $1 ; shift ; done

CWD=`pwd`

############################ SANITY CHECKS ########################

# Check for a flash.par.
[ -f $FLASH_PAR ] || {

  echo "There is no file $FLASH_PAR to use as a flash.par template.  Please"
  echo "supply one, or point to one on the command line by setting FLASH_PAR."
  exit 1
}

# Check for the source tree
[ -d $FLASH_PATH ] || {

  echo "$FLASH_PATH is not a directory or a link to a directory.  As a"
  echo "consequence there's no Flash source tree to work with.  Please"
  echo "correct this."
}

# Check for the build directory
[ -d $FLASH_PATH/$FLASH_BUILD_DIR ] || {

  echo "No Flash build directory was found at $FLASH_PATH/$FLASH_BUILD_DIR."
  echo "That's where we get the executable and its files, you know."
  exit 1
}

# Let's find out how many properly-formed run directories with this
# stem already exist.
stlen_plus_one=`echo $RUNDIR_STEM|wc|awk '{print $3}'`
last_rundir=`ls | grep "$RUNDIR_STEM[0-9][0-9][0-9]$" | tail -1`

[ "$last_rundir" ] && {

  seqnum_prev=`echo $last_rundir | cut -c ${stlen_plus_one}-`

} || {

  seqnum_prev=0
}

# Smartass alert
[ $seqnum_prev -eq 999 ] && {

  echo "Apparently there is a run directory called ${RUNDIR_STEM}_999. There"
  echo "is nowhere for a dumb script like this to go to increment the run"
  echo "directory sequence number.  You supply the required intelligence."
  exit 1
}

# Check than any requests for restarts are sane.
if [ \( "$RESTART_CHK" -a ! "$RESTART_LOG" \) -o \
       \( ! "$RESTART_CHK" -a "$RESTART_LOG" \) ]; then

 echo "For a restart, you must specify both a checkpoint file (by setting RESTART_CHK)"
 echo "and a log file (by specifying RESTART_LOG)."
 exit 1

fi

if [ "$RESTART_CHK"  ]; then

 ls $RESTART_CHK >& /dev/null || {  # This should be '[ -f $RESTART_CHK ]', but
                                    # AIX is brain-dead, and test doesn't work
                                    # with large files.

    echo "$RESTART_CHK not found, regrets."
    exit 1
  }
fi

if [ "$RESTART_LOG"  ]; then

  [ -f $RESTART_LOG ] || {

    echo "$RESTART_LOG not found, regrets."
    exit 1
  }
fi

if [ "$EXISTING_RUNDIR" ]; then

  rd=`printf "%s%3.3d" $RUNDIR_STEM $EXISTING_RUNDIR`
  [ -d $rd ] || {

    echo "Alleged existing run directory $rd not found, regrets."
    exit 1
  }
fi

################### BASIC DIRECTORY CREATION/POPULATION ########################

if [ "$EXISTING_RUNDIR" ]; then

  new_rundir=$rd

else

# Make the new run directory
  seqnum=`echo $seqnum_prev | awk '{print $1 + 1}'`
  new_rundir=`printf "%s%3.3d" $RUNDIR_STEM $seqnum`
  echo "Creating new run directory $new_rundir..."
  mkdir "$new_rundir" || {

    echo "Directory creation failed.  That was unexpected."
    exit 1
  }

fi

# It needs an executable, and any required data files.
[ -x $FLASH_PATH/$FLASH_BUILD_DIR/flash4 ] || {

  echo "No flash4 executable found in $FLASH_PATH/$FLASH_BUILD_DIR.  Did you"
  echo "build one?"
  [ ! "$EXISTING_RUNDIR" ] && rmdir $new_rundir
  exit 1
}

echo "Copying flash4 and any required datafiles to $new_rundir..."
cp $FLASH_PATH/$FLASH_BUILD_DIR/flash4 $new_rundir
[ -f $FLASH_PATH/$FLASH_BUILD_DIR/setup_datafiles ] && {

  for path in `cat $FLASH_PATH/$FLASH_BUILD_DIR/setup_datafiles`; do
    datafile=`basename $path`
    cp $FLASH_PATH/$FLASH_BUILD_DIR/$datafile $new_rundir
  done
}

[ -f $BATCH_SCRIPT ] && {

  cp $BATCH_SCRIPT $new_rundir
  echo "Copied $BATCH_SCRIPT to $new_rundir"

} || {

  echo "Batch submission script not found, bring your own if you need one."
}

touch $new_rundir/OPERATOR_NOTES

################### SOURCE VERSION INFO ########################

# Let's get version information
echo "Saving version information.  This might take a minute (it's a big tree...)"
( cd $FLASH_PATH/source ; svn info ) > $new_rundir/FLASH_SVN_INFO || {

  echo "Subversion apparently doubts that $FLASH_PATH is a genuine, uncorrupted,"
  echo "svn-versioned source directory, produced by a 'svn checkout'.  It is generally"
  echo "right about that sort of thing.  So please supply a source directory that meets"
  echo "its expectations."
  [ ! "$EXISTING_RUNDIR" ] && rm -rf $new_rundir
  exit 1
}

( cd $FLASH_PATH/source ; echo "" ; svn status --verbose ) >> $new_rundir/FLASH_SVN_INFO
changed=`tail -n +12 $new_rundir/FLASH_SVN_INFO | grep '^[ADMR]' | awk '{print $NF}'`
( cd $FLASH_PATH/source ; svn diff $changed ) > $new_rundir/FLASH_SVN_DIFF

[ -s $new_rundir/FLASH_SVN_DIFF ] && {

  echo "Note:  there are local modifications of the source tree.  Details of"
  echo "the modifications have been saved in $new_rundir/FLASH_SVN_DIFF."

} || {

  rm $new_rundir/FLASH_SVN_DIFF
}
echo "Version information locked in."


################### FUN WITH flash.par ########################

# Munge the base name
echo "Tweaking flash.par..."
basebasename=`grep '^basenm' $FLASH_PAR | cut -d= -f2 | tr -d \"`
[ "$basebasename" ] &&
  basename=`printf "%s_%s" $new_rundir $basebasename` ||
  basename=$new_rundir

if [ ! "$RESTART_CHK" ] ; then    # New start



 # (make sure basenm is there)
    cp $FLASH_PAR flash.par.tmp
    grep '^basenm' flash.par.tmp >& /dev/null ||
     echo 'basenm' >> flash.par.tmp 

 cat flash.par.tmp | 

 #       -- munge the base name...
    sed -e "/^basenm/c\\
basenm		= \"$basename\""  |

 #       -- make sure restart = .false. ...

    sed -e "/^restart/c\\
restart		= .false." |

 #       -- comment out any setting of plotFileNumber and ParticleFileNumber

    sed -e "s/^plotFileNumber/# plotFileNumber/" \
        -e "s/^ParticleFileNumber/# ParticleFileNumber/"  > $new_rundir/flash.par

    rm flash.par.tmp


else      # Restart, have checkpoint file and log file

 

  # Get index of checkpoint file
  chk_ind_pos=`echo $RESTART_CHK | awk -F _ '{print NF}'`
  chk_ind=`echo $RESTART_CHK | cut -d_ -f $chk_ind_pos`
  [ "$chk_ind" ] && [ ! "`echo $chk_ind | tr -d [:digit:]`" ] || {

    echo "The alleged checkpoint file $RESTART_CHK has a suspicious filename -- it doesn't"
    echo "end in _XXX, where XXX are digits.  This is a requirement, sorry."
    [ ! "$EXISTING_RUNDIR" ] && rm -rf $new_rundir
    exit 1
  }

  chk_name=${basename}${CHK_SUFF}$chk_ind
  [ -e $new_rundir/$chk_name ] && rm $new_rundir/$chk_name
  if [ "`echo $RESTART_CHK | cut -c1`" = "/" ]; then

    ln -s $RESTART_CHK $new_rundir/$chk_name  # $RESTART_CHK is a full pathname

  else

    ln -s ../$RESTART_CHK $new_rundir/$chk_name # $RESTART_CHK is a relative pathname

  fi

  # Get indices of plot and particle files corresponding to this restart.
  # Got to make sure this log file knows about this checkpoint file first:
  chk_base=`basename $RESTART_CHK`
  ( grep $chk_base $RESTART_LOG > /dev/null ) || {
    echo "The log file $RESTART_LOG has never heard of the checkpoint file $RESTART_CHK."
    echo "That means the plot and particle index numbers cannot be adjusted properly, alas."
    [ ! "$EXISTING_RUNDIR" ] && rm -rf $new_rundir
    exit 1
  }

  # plotfile index
  plt_name=`cat $RESTART_LOG |
            awk '/\[IO_writePlotfile\] close: type=plotfile/ {plt=$8}
                 /\[IO_writeCheckpoint\] close: type=checkpoint name='$chk_base'/ {print plt ; exit}'`

  plt_ind_pos=`echo $plt_name | awk -F _ '{print NF}'`
  plt_ind=`echo $plt_name | cut -d_ -f $plt_ind_pos`
  echo "Last Plotfile index: $plt_ind."
  plt_ind_new=`echo $plt_ind | awk '{print $1 + 1}'`
  echo "Will start plots from index $plt_ind_new."

  # particle file index
  part_name=`cat $RESTART_LOG |
            awk '/\[IO_writeParticles\] close: type=particles/ {part=$8}
                 /\[IO_writeCheckpoint\] close: type=checkpoint name='$chk_base'/ {print part ; exit}'`

  if [ "$part_name" ]; then
  
    part_ind_pos=`echo $part_name | awk -F _ '{print NF}'`
    part_ind=`echo $part_name | cut -d_ -f $part_ind_pos`
    echo "Last particle file index: $part_ind."
    part_ind_new=`echo $part_ind | awk '{print $1 + 1}'`
    echo "Will start particle files from index $part_ind_new."

  else
  
    echo "No particles were written according to the logfile, so we won't"
    echo "adjust the particle file number."

  fi

# Hack to make sure the required keywords exist in the flash.par file
  cp $FLASH_PAR flash.par.tmp
  grep '^basenm' flash.par.tmp >& /dev/null ||
     echo "basenm" >> flash.par.tmp

  grep '^restart' flash.par.tmp >& /dev/null ||
     echo "restart" >> flash.par.tmp

  grep '^checkpointFileNumber' flash.par.tmp >& /dev/null ||
    echo "checkpointFileNumber" >> flash.par.tmp

  grep '^plotFileNumber' flash.par.tmp >& /dev/null ||
    echo "plotFileNumber" >> flash.par.tmp

  cat flash.par.tmp |

 #       -- munge the base name...

    sed -e "/^basenm/c\\
basenm		= \"$basename\""  |

 #       -- make sure restart = .true. ...

    sed -e "/^restart/c\\
restart		= .true." |

 #       -- set the checkpoint file number...

    sed -e "/^checkpointFileNumber/c\\
checkpointFileNumber	= $chk_ind"  |

 #       -- set the plot file number...
 
    sed -e "/^plotFileNumber/c\\
plotFileNumber		= $plt_ind_new"  > $new_rundir/flash.par

  rm flash.par.tmp

 # --- and the particle file number, if required.

  if [ "$part_name" ]; then

    mv $new_rundir/flash.par flash.par.tmp
    grep '^ParticleFileNumber' flash.par.tmp >& /dev/null ||
      echo "ParticleFileNumber" >> flash.par.tmp
    cat flash.par.tmp |
      sed -e "/^ParticleFileNumber/c\\
ParticleFileNumber	= $part_ind_new"  > $new_rundir/flash.par

    rm flash.par.tmp

  fi

fi

