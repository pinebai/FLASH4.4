#!/bin/bash -x
HPCTOOLKIT_PATH='/home/projects/hpctoolkit/pkgs/hpctoolkit'
HPCTOOLKIT_MEASUREMENTS='hpctoolkit-flash4-measurements-484550'
FLASH_OBJ='/intrepid-fs0/users/cdaley/persistent/2012/multithreaded/bgq/rtflame_maxblocks_120_tbl_hpctoolkit_par'

qsub -A TurbNuclComb_esp -q prod-devel --env BG_COREDUMPONERROR=1:MPIRUN_LABEL=1 -t 1:00:00 -n 512 --mode smp \
    ${HPCTOOLKIT_PATH}/bin/hpcprof-mpi -S flash4.hpcstruct -I "${FLASH_OBJ}/*" ${HPCTOOLKIT_MEASUREMENTS}
