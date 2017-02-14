/* Each file that wants to use global variables must include io_ncmpi_use_globals.h
   header file.  The only file that does not include io_ncmpi_use_globals.h is
   io_ncmpi_set_globals.h where the global variables are defined. */

#ifndef IO_SET_GRID_MPI_TYPES_H
#define IO_SET_GRID_MPI_TYPES_H

#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include "constants.h"
#include "io_flash.h"
#include "mangle_names.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

void FTOC(io_init_grid_mpi_types)(const int * const pMyPE,
				  const int * const pGridDataStruct,
				  const int blockOuterSize[],
				  const int blockInnerSize[],
				  const int blockInnerOffset[],
				  const int * const pNumGridVar,
				  const int plotVarArr[],
				  const int * const pNumPlotVar);

void FTOC(io_free_grid_mpi_types)(void);


int io_resize_mpi_type(MPI_Datatype oldtype,
		       MPI_Aint lb,
		       MPI_Aint extent,
		       MPI_Datatype *pNewtype);

#endif
