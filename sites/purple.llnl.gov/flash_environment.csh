# Environment settings for University of Chicago Flash Center operations

umask 07

setenv FLASH_GROUP_ID flashers

setenv FLASH_PROJECT_HOME ~carlo/project

setenv FLASH_SOURCE_TREE $FLASH_PROJECT_HOME/Flash3

setenv FLASH_SITE_NAME purple.llnl.gov

set path =  ($path $FLASH_PROJECT_HOME/bin)

setenv MP_COREFILE_FORMAT 1   # Generate lightweight core file format,
                              # to file rather than to stdout

source /usr/global/tools/dotkit/init.csh
use svn
