#!/bin/bash -x

# Notes
#
# 1). Sub-block jobs < 128 nodes do not work properly yet.
#     It is best to avoid the trouble and just use 128 nodes.
# 2). I get failures in malloc when using BG_STACKGUARDENABLE=1
#     I think this is a system error and not the stack and heap
#     running into each other.

PROG=flash4
NODES=128
RANKS_PER_NODE=16
NPROCS=$((NODES*RANKS_PER_NODE)) 

#BG_THREADLAYOUT:  1 - default next core first; 2 - my core first
qsub -A TurbNuclComb_esp -n ${NODES} --mode c${RANKS_PER_NODE} -t 0:30:00 \
 --env BG_THREADLAYOUT=2:BG_SHAREDMEMSIZE=32:BG_COREDUMPONERROR=1:OMP_NUM_THREADS=4:OMP_STACKSIZE=16M:L1P_POLICY=std:PAMID_VERBOSE=1 \
 $PROG
