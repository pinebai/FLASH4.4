#!/bin/bash -x

HPCTOOLKIT_PATH='/home/projects/hpctoolkit/pkgs/hpctoolkit'
HPCTOOLKIT_MEASUREMENTS='hpctoolkit-flash4-measurements-492281'
FLASH_OBJ='/intrepid-fs0/users/cdaley/persistent/2012/multithreaded/bgq/rtflame_maxblocks_120_hpctoolkit'

#Create flash4.hpcstruct
#${HPCTOOLKIT_PATH}/bin/hpcstruct --loop-fwd-subst=no ./flash4

#I have found that hpcprof-mpi can consume more memory than is available
#in vn mode and so I use smp mode instead.
qsub -A TurbNuclComb_esp -q prod-devel --env BG_COREDUMPONERROR=1:MPIRUN_LABEL=1 -t 1:00:00 -n 512 --mode smp \
    ${HPCTOOLKIT_PATH}/bin/hpcprof-mpi -S flash4.hpcstruct -I "${FLASH_OBJ}/*" ${HPCTOOLKIT_MEASUREMENTS}
