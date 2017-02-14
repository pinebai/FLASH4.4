#!/bin/sh
#This script tells hpctoolkit to use delayed sampling.  This will only work
#if you have built flash with monitors/Profiler/ProfilerMain/hpctoolkit unit.
#We make use of delayed sampling to measure events in the evolution section
#of FLASH only.
qsub -A TurbNuclComb_esp -q prod-devel --env BG_MAPPING=TXYZ:BG_COREDUMPONERROR=1:MPIRUN_LABEL=1:HPCRUN_EVENT_LIST="WALLCLOCK@5000":HPCRUN_DELAY_SAMPLING=1 -t 1:00:00 -n 512 --mode smp `pwd`/flash4
