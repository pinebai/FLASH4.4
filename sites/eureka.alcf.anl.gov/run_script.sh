#!/bin/bash
NODES=1
RANKS_PER_NODE=8
NPROCS=$((NODES*RANKS_PER_NODE))
DEFAULT_ENV="mx"
BINARY="./flash4"


#Bash parameter substitution - use the default environment if no argument.
MY_ENV=${1:-$DEFAULT_ENV}
echo -e "Using software environment ${MY_ENV}\n\
Running ${BINARY} with ${NPROCS} MPI ranks \
(${NODES} nodes and ${RANKS_PER_NODE} MPI ranks per node)"
#We want resoft function available in our environment.
source /etc/profile.d/00softenv.sh
#We use a custom .softenvrc file which makes use of MY_ENV variable.
resoft
echo -e "PATH is $PATH\nLD_LIBRARY_PATH is $LD_LIBRARY_PATH\n"


if [ "$MY_ENV" = "mx" ]; then
    # For mpd process manager.
    # Includes +mpich2-mx-1.0.7..2 and +mpich2-mx-1.2.1p1..8-intel.
    mpdboot -n $NODES -f $COBALT_NODEFILE
    mpiexec -n $NPROCS $BINARY
elif [ "$MY_ENV" = "nomx" ]; then
    # For hydra process manager.
    # Includes +mpich2-1.3.1-intel.
    mpiexec -n $NPROCS -f $COBALT_NODEFILE $BINARY
else
    echo "Environment not recognised!"
fi
