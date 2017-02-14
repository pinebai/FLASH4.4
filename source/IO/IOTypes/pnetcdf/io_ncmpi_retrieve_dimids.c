#include <assert.h>
#include <string.h>
#include <stdio.h>

#include <mpi.h>
#include <pnetcdf.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "mangle_names.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

/* Retrieves an array of dimension IDs that are required to
   specify the size of the dataset in the file */

void FTOC(io_ncmpi_retrieve_dimids)(const int * const pMyPE,
				    const int * const pFileID,
				    const int * const pFileFmt,
				    const int * const pGridStruct,
				    int dimIDs[IO_MESH_DIMS])
{
  /* arrays of IO_MESH_DIMS pointers to const char */
  const char * unk[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb","dim_nUnkVar"};
  const char * facex[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb_long","dim_nFacexVar"};
  const char * facey[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb_long","dim_nxb","dim_nFaceyVar"};
  const char * facez[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb_long","dim_nyb","dim_nxb","dim_nFacezVar"};
  const char * scratch[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb","dim_nScratchGridVar"};

  const char * unk_single_var[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb",NULL};
  const char * facex_single_var[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb_long",NULL};
  const char * facey_single_var[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb_long","dim_nxb",NULL};
  const char * facez_single_var[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb_long","dim_nyb","dim_nxb",NULL};
  const char * scratch_single_var[IO_MESH_DIMS] =
    {"dim_tot_blocks","dim_nzb","dim_nyb","dim_nxb",NULL};

  /* pointer to an array of IO_MESH_DIMS pointers to const char... YUK! */
  const char * (*pDataStructure)[IO_MESH_DIMS] = NULL;
  const char * pDimName = NULL;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int fileFmt = *pFileFmt;
  const int gridStruct = *pGridStruct;
  int i, err;
  MPI_Offset len;


  switch(gridStruct)
  {
  case CENTER:
    if (fileFmt == 9) {
      pDataStructure = &unk_single_var;
    } else {
      pDataStructure = &unk;
    }
    break;
  case SCRATCH:
    if (fileFmt == 9) {
      pDataStructure = &scratch_single_var;
    } else {
      pDataStructure = &scratch;
    }
    break;
  case FACEX:
    if (fileFmt == 9) {
      pDataStructure = &facex_single_var;
    } else {
      pDataStructure = &facex;
    }
    break;
  case FACEY:
    if (fileFmt == 9) {
      pDataStructure = &facey_single_var;
    } else {
      pDataStructure = &facey;
    }
    break;
  case FACEZ:
    if (fileFmt == 9) {
      pDataStructure = &facez_single_var;
    } else {
      pDataStructure = &facez;
    }
    break;
  default:
    Driver_abortFlashC("[io_ncmpi_retrieve_dimids] Unknown Grid structure");
  }


  /* Retrieve the dimension IDs representing global size in each dimension */
  for (i=0; i<IO_MESH_DIMS; ++i) {
    pDimName = (*pDataStructure)[i];

    if (debugIO) {
      if (myPE == MASTER_PE) {
	printf(" [io_ncmpi_retrieve_dimids]: Getting dimID[%d] which "
	       "corresponds to name %s\n", i, pDimName);
      }
    }

    if (pDimName) {
      err = ncmpi_inq_dimid(fileID, pDimName, &dimIDs[i]);
      assert(err == NC_NOERR);

      if (debugIO) {
	err = ncmpi_inq_dimlen(fileID, dimIDs[i], &len);
	assert(err == NC_NOERR);
	if (myPE == MASTER_PE) {
	  printf(" [io_ncmpi_retrieve_dimids]: dimID[%d]=%lld\n", i, len);
	}
      }
    } else {
      dimIDs[i] = -1;
    }
  }
}
