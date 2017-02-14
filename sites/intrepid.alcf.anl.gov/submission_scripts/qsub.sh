#!/bin/sh
#Important.  We must use XLSMPOPTS=stack=N instead of OMP_STACKSIZE=N.  This is
#because OMP_STACKSIZE is ignored by the OpenMP implementation on BG/P since
#the OpenMP version on BG/P is OpenMP-2.5 (200505).  The OMP_STACKSIZE
#environmental variable was introduced in OpenMP-3.0.
qsub -A TurbNuclComb_esp -q prod-devel --env BG_MAPPING=TXYZ:BG_COREDUMPONERROR=1:MPIRUN_LABEL=1:OMP_NUM_THREADS=4:XLSMPOPTS=stack=16000000 -t 1:00:00 -n 256 --mode smp `pwd`/flash4
