#!/bin/bash -x
#
# Asssumptions
#  * the qsub line requests a script mode, e.g.
#    qsub ... --mode script cobalt_script.sh
#  * the data files and flash4 binary are in the current working directory
#  * the thread experiment directories do not exist
#
# Notes
#  * COBALT_JOBSIZE shows the number of nodes requested in the qsub line
#  * LOCARGS allows us to use this script for sub-block jobs (< 128 nodes) too

PROG=flash4
NODES=${COBALT_JOBSIZE}
RANKS_PER_NODE=16
RUNFILES="${PROG}
flash.par
amr_runtime_parameters
helm_table.dat
coldwd_mchandra_damped.dat
nse_dens_emq_table.txt
nse_pres_hmq_table.txt
SpeciesList.txt"

NPROCS=$((NODES*RANKS_PER_NODE))
LOCARGS="--block $COBALT_PARTNAME ${COBALT_CORNER:+--corner} $COBALT_CORNER ${COBALT_SHAPE:+--shape} $COBALT_SHAPE"

for i in 1 2 3 4; do

    echo "Running ${PROG} with ${NPROCS} MPI processes and ${i} OpenMP threads"
    OUTPUT=N${NODES}_R${RANKS_PER_NODE}_T${i}_${PROG}_${NPROCS}_ranks.log
    mkdir ${i}
    cd ${i}
    for j in ${RUNFILES}; do
	ln -s ../${j}
    done

    #BG_THREADLAYOUT:  1 - default next core first; 2 - my core first
    runjob \
	--ranks-per-node ${RANKS_PER_NODE} \
	--np ${NPROCS} \
	--cwd ${PWD} \
	${LOCARGS} \
	--envs BG_THREADLAYOUT=2 \
	--envs BG_SHAREDMEMSIZE=32 \
	--envs BG_COREDUMPONERROR=1 \
	--envs OMP_NUM_THREADS=${i} \
	--envs OMP_STACKSIZE=16M \
	--envs L1P_POLICY=std \
	--envs PAMID_VERBOSE=1 \
	: ${PROG} > $OUTPUT 2>&1
    cd ..
done
