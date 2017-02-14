#include "io_use_grid_mpi_types.h"
#include <stdio.h>

/* Returns a pointer to the global MPI type describing the
   particular file type / grid type combination */
MPI_Datatype * io_get_grid_mpi_type(const int fileFmt,
				    const int fileType,
				    const int gridDataStruct)
{
  MPI_Datatype *pMPI_DataType = NULL;

  assert(fileFmt <= 10);
  assert(gridDataStruct >= 1 && gridDataStruct <= 5);

  if (fileFmt <= 9) {
    pMPI_DataType = &G_mesh_single_var_type[gridDataStruct-1];
  } else if (fileFmt == 10) {

    switch(fileType)
    {
    case CHECKPOINTFILE:
      pMPI_DataType = &G_mesh_all_var_type[gridDataStruct-1];
      break;
    case PLOTFILE:
      pMPI_DataType = &G_mesh_subset_var_type[gridDataStruct-1];
      break;
    default:
      Driver_abortFlashC("Unexpected filetype");
    }
  }

  assert(*pMPI_DataType != MPI_DATATYPE_NULL);
  return pMPI_DataType;
}
