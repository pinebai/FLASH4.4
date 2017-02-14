# Environment settings for University of Chicago Flash Center operations

umask 07

export FLASH_GROUP_ID=TurbNuclComb_esp
export FLASH_SITE_NAME=vesta.alcf.anl.gov

export FLASH_PROJECT_HOME=/gpfs/vesta-fs0/projects/TurbNuclComb_esp/${USER:-${LOGNAME}}/2013
export FLASH_SOURCE_TREE=/home/${USER:-${LOGNAME}}/flash/laser_performance

# Add coreprocessor.pl and bgq_stack to our PATH
PATH=${PATH}:/soft/debuggers/scripts/bin

# Add HDF5 tools such as h5dump and h5diff to our PATH.
HDF5_PATH=/gpfs/vesta-fs0/projects/TurbNuclComb_esp/cdaley/software/current/hdf5/1.8.10/fen/xl
PATH=${PATH}:${HDF5_PATH}/bin
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HDF5_PATH}/lib

# Add sfocu to our PATH.
#  * The binary named sfocu compares grid data only.
#  * The binary named sfocu.particles compares grid and particles data.
#    - comparing N particles can be extremely slow becuase there is an
#      N^{2} particle comparison algorithm.
PATH=${PATH}:/gpfs/vesta-fs0/projects/TurbNuclComb_esp/cdaley/software/current/sfocu/fen/xl

export PATH
export LD_LIBRARY_PATH

# Your ~/.soft file should contain
# +mpiwrapper-xl
# @default
