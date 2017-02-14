#!/bin/bash
qsub -q prod-devel -A TurbNuclComb_esp -t 60 -n 256 --env BG_MAPPING=TXYZ:BG_COREDUMPONERROR=1:MPIRUN_LABEL=1:BG_MAXALIGNEXP=0:XLSMPOPTS=stack=16000000 --mode script cobalt_job.sh
