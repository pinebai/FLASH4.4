#ifndef IO_USE_GRID_MPI_TYPES_H
#define IO_USE_GRID_MPI_TYPES_H

#include <mpi.h>
#include <assert.h>
#include "constants.h"
#include "io_flash.h"

/* MPI variable that stores datatype describing UNK block interior */

extern MPI_Datatype G_mesh_all_var_type[NUM_MESH_TYPES],
  G_mesh_subset_var_type[NUM_MESH_TYPES],
  G_mesh_single_var_type[NUM_MESH_TYPES];

MPI_Datatype * io_get_grid_mpi_type(const int fileFmt,
				    const int fileType,
				    const int gridDataStruct);
#endif
