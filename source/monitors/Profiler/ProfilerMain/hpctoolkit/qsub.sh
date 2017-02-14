#!/bin/bash -x
#Note HPCRUN_DELAY_SAMPLING=1 which tell HPCToolkit that we are enabling
#sampling for a subset of our application.
qsub -A TurbNuclComb_esp -q prod-devel --env BG_MAPPING=TXYZ:BG_COREDUMPONERROR=1:MPIRUN_LABEL=1:OMP_NUM_THREADS=4:XLSMPOPTS=stack=16000000:HPCRUN_EVENT_LIST="WALLCLOCK@5000":HPCRUN_DELAY_SAMPLING=1 -t 1:00:00 -n 256 --mode smp `pwd`/flash4
