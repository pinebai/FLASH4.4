# Environment settings for University of Chicago Flash Center operations

umask 07

export FLASH_GROUP_ID=incite13
export MP_COREFILE_FORMAT=1   # Generate lightweight core file format,
                              # to file rather than to stdout
export FLASH_SITE_NAME=seaborg.nersc.gov

export FLASH_PROJECT_HOME=/project/projectdirs/incite13/project/seaborg
export FLASH_SOURCE_TREE=$FLASH_PROJECT_HOME/Flash3_wd_def

export PATH=${PATH}:$FLASH_PROJECT_HOME/bin

module load python
