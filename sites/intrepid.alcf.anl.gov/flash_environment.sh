# Environment settings for University of Chicago Flash Center operations

umask 07

export FLASH_GROUP_ID=flash
export MP_COREFILE_FORMAT=1   # Generate lightweight core file format, -- from Seaborg, may be useless here
                              # to file rather than to stdout.         -- from Seaborg, may be useless here 
#export BG_COREDUMPONEXIT=1   # forces all processors to dump core on exit.
export FLASH_SITE_NAME=intrepid.alcf.anl.gov

export FLASH_PROJECT_HOME=/gpfs1/${USER:-${LOGNAME}}/2009
export FLASH_SOURCE_TREE=$FLASH_PROJECT_HOME/src/project-WD_paramStudy

pathmunge () {
        if ! echo $PATH | /bin/egrep -q "(^|:)$1($|:)" ; then
           if [ "$2" = "after" ] ; then
              PATH=$PATH:$1
           else
              PATH=$1:$PATH
           fi
        fi
}

pathmunge /soft/apps/hdf5-1.8.0-fen/bin after
pathmunge $FLASH_PROJECT_HOME/bin after

unset pathmunge

