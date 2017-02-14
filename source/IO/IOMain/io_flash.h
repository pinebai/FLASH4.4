#if 0
This header file must be kept usable by C and Fortran code.
#endif

#ifndef IO_FLASH_H
#define IO_FLASH_H

#include "Flash.h"

#ifdef FLASH_GRID_UG
#ifdef FIXEDBLOCKSIZE
#define IO_FLASH_UG
#else
#define IO_FLASH_NOFBS_UG
#endif
#endif

#define IO_FLASH_INT 1
#define IO_FLASH_DOUBLE 2
#define IO_FLASH_FLOAT 3
#define IO_FLASH_CHAR 4
#define IO_FLASH_STRING 5

#define IO_READ_XFER 1
#define IO_WRITE_XFER 2
#define IO_READ_XFER_MASTER_PE 3
#define IO_WRITE_XFER_MASTER_PE 4

#define IO_FILE_HDF5 1
#define IO_FILE_PNETCDF 2

#define IO_MAX_DIMS 5
#define IO_MESH_DIMS 5
#define NUM_MESH_TYPES 5

#define LEGACY_LABEL_LEN 4

#ifdef __FLASH_HEADER_FILE__
# if (NUNK_VARS > NFACE_VARS)
#  if (NUNK_VARS > NSCRATCH_GRID_VARS)
#   define MAX_VAL (NUNK_VARS)
#  else
#   define MAX_VAL (NSCRATCH_GRID_VARS)
#  endif
# else
#  if (NFACE_VARS > NSCRATCH_GRID_VARS)
#   define MAX_VAL (NFACE_VARS)
#  else
#   define MAX_VAL (NSCRATCH_GRID_VARS)
#  endif
# endif
#else
# define MAX_VAL 100
#endif

#if 0
The code below is just to make type based I/O compatible with
mesh replication.  Many arrays are statically sized in type based
I/O because I assumed that we would always know the number of
variables in a simulation at compile time.  This stopped being
the case when mesh replication was introduced.  The hack below
gives us more space to play with when using simulations with
mesh replication (we check for possible overflow in io_typeInit).
Type based I/O needs to be refactored in more ways than making
a few arrays dynamic.
#endif
#if NONREP_COUNT > 0
#define MAX_MESH_VAR (MAX_VAL * 3)
#else
#define MAX_MESH_VAR (MAX_VAL)
#endif


#define CHECKPOINTFILE 888
#define PLOTFILE 889


#if 0
Print debugging information - Need to move DEBUG_IO defintion.
#define DEBUG_IO
#endif

#ifdef PROFILE_IO
#define IO_TIMERS_START(s) call Timers_start(s)
#define IO_TIMERS_STOP(s) call Timers_stop(s)
#else
#define IO_TIMERS_START(s)
#define IO_TIMERS_STOP(s)
#endif

#define IO_CHECK_XFER(rtn,dataset) \
if (rtn/=0) call Driver_abortFlash \
('[' // FILE_AT_LINE // '] Transfer error: ' // dataset)

#endif
