# Environment settings for University of Chicago Flash Center operations

umask 07

export FLASH_GROUP_ID=flashers

export FLASH_PROJECT_HOME=~carlo/project

export FLASH_SOURCE_TREE=$FLASH_PROJECT_HOME/Flash3

export FLASH_SITE_NAME=purple.llnl.gov

export PATH=${PATH}:$FLASH_PROJECT_HOME/bin

export MP_COREFILE_FORMAT=1   # Generate lightweight core file format,
                              # to file rather than to stdout

. /usr/global/tools/dotkit/init.sh
use svn
