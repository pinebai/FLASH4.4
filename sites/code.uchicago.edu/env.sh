#!/bin/bash
#
# Usage:
# > . env.sh gnu
# > . env.sh intel
# > . env.sh sun
# (Choose one)
#
# This script modifies environmental variables for gnu, intel, sun software stacks.

thisfile="$BASH_SOURCE"

if [ "$1" != "gnu" -a "$1" != "intel" -a "$1" != "intel-static" -a "$1" != "sun" ]; then
  echo "ERROR: Unknown switch '$1'. Accepted values: gnu, intel, sun, intel-static"
  return 1;
fi
if [ "$1" == "intel-static" ]; then
    STACK="${1%-static}"	# strip off the -static
else
    STACK=${1}
fi
staticflag=""		# reset flag
if [ "$1" != "${STACK}" ]; then
    staticflag="y"		# -static was stripped off
fi

# Source the definition of pathmungeany()
# NOTE: CHANGE THIS IF THIS FILE (sites/<SOMESITE>/env.sh) IS COPIED FROM ITS
# LOCATION IN THE SOURCE TREE.

source $(dirname "$thisfile")/../../tools/scripts/fcndefs.bash


if [ "${STACK}" == "gnu" ]; then
  # So that gfortran does not buffer stdout:
  # This may slow down the application slightly, but without this setting
  # we only get FLASH stdout every few hundred time steps.
  export GFORTRAN_UNBUFFERED_ALL='y'
elif [ "${STACK}" == "intel" ]; then
  # So that ifort compiled applications dump core:
  export decfort_dump_flag='y'
elif [ "${STACK}" == "sun" ]; then
  #Add the oracle studio compiler suite.
  pathmungeany /opt/solstudio12.3/bin first
  #Enable OpenMP runtime checking.
  export SUNW_MP_WARN=TRUE
fi


# Create core dumps:
ulimit -c unlimited

# So that electric fence does not abort in MPI_init:
export EF_ALLOW_MALLOC_0=1


# Usable by gnu and intel (even though compiled with gnu):
# ----------------------------------------------------------------------

VALGRIND=/opt/valgrind
pathmungeany ${VALGRIND}/bin

# So we can produce a valgrind log file for each mpi process.
# Usage:
# > mpiexec -n 4 $valgrind_log ./flash3
export valgrind_log='valgrind --tool=memcheck --log-file=valgrind.log.%p --time-stamp=yes --track-origins=yes --run-libc-freeres=yes'

# Add fancy archer branch of GDB to the path:
#NOT YET AVAILABLE
#export PATH="/usr/local/gdb-archer/bin:$PATH"

# 6.3 Building and Installing HPCToolkit
# Q: Do I need to compile HPCToolkit with any special options for MPI
# support?
# A: No, HPCToolkit is designed to work with multiple MPI
# implementations at the same time. That is, you don’t need to provide
# an mpi.h include path, and you don’t need to compile multiple versions
# of HPCToolkit, one for each MPI implementation.  The
# technically-minded reader will note that each MPI implementation uses
# a different value for MPI_COMM_WORLD and may wonder how this is
# possible. hpcrun (actually libmonitor) waits for the application to
# call MPI_Comm_rank() and uses the same communicator value that the
# application uses. This is why we need the application to call
# MPI_Comm_rank() with communicator MPI_COMM_WORLD.
HPCTOOLKIT=/home/cdaley/software/hpctoolkit/r3542
pathmungeany ${HPCTOOLKIT}/bin

# Example usage of hpctoolkit to profile FLASH (no instrumentation necessary):
# > hpcrun ./flash3
# > hpcstruct flash3
# > hpcprof -S flash3.hpcstruct -I /home/cdaley/flash/trunk/jeans hpctoolkit-flash3-measurements
# (The -I is the object directory containing the source files of
# the profiled flash3 binary)

# View results using hpcviewer
# > hpcviewer hpctoolkit-flash3-database


# NOTE: Only usable in parallel by gnu's MPI.
# I've not compiled a parallel intel version of visit yet.
# > visit -np 8
if [ "${STACK}" == "gnu" ]; then
  VISIT=/opt/visit/2.3.2
  pathmungeany ${VISIT}/bin after
fi


# Select the correct MPI, HDF5 and Pnetcdf for gnu or intel:
# ----------------------------------------------------------------------
MPI=/opt/mpich2/${STACK}/1.4.1p1

if [ -z  "$staticflag" ]; then
  pathmungeany ${MPI}/lib first LD_LIBRARY_PATH
else
  pathmungeany ${MPI}/lib first LIBRARY_PATH # Is LIBRARY_PATH actually used??
fi

HDF5=/opt/hdf5/${STACK}/1.8.7

PNET=/opt/netcdf/${STACK}/1.2.0

pathmungeany ${MPI}/bin first
pathmungeany ${HDF5}/bin
pathmungeany ${PNET}/bin

pathmungeany ${PNET}/man:${MPI}/share/man first MANPATH


export PATH
export LD_LIBRARY_PATH
