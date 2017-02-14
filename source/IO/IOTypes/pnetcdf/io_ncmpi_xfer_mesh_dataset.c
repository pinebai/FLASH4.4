#include "io_ncmpi_xfer_mesh_dataset.h"

void io_ncmpi_xfer_mesh_dataset(const int myPE,
				const int fileID,
				const int xferType,
				const int fileType,
				const int fileFmt,
				const int gridDataStruct,
				const int numDataDims,
				const int nonBlocking,
				const int globalOffset[],
				const int localSubSize[],
				const char datasetName[],
				double * pData)
{
  FLASH_IO_Request *pNCMPI_Request = NULL;
  MPI_Datatype *pMPI_DataType = NULL;
  MPI_Offset startDisk[IO_MESH_DIMS], count[IO_MESH_DIMS];
  int i, err, varID;
  char *fileTypeText = NULL, *xferTypeText = NULL;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  assert(myPE >= 0);
  assert(xferType == IO_WRITE_XFER || xferType == IO_READ_XFER);
  assert(fileType == CHECKPOINTFILE || fileType == PLOTFILE);
  assert(gridDataStruct == CENTER || gridDataStruct == FACEX ||
	    gridDataStruct == FACEY || gridDataStruct == FACEZ ||
	    gridDataStruct == SCRATCH);
  assert(numDataDims <= IO_MESH_DIMS);

  if (debugIO) {
    if (fileType == CHECKPOINTFILE) fileTypeText = "checkpoint file";
    if (fileType == PLOTFILE) fileTypeText = "plot file";
    if (xferType == IO_WRITE_XFER) xferTypeText = "write";
    if (xferType == IO_READ_XFER) xferTypeText = "read";
    if (myPE == MASTER_PE) {
      printf(" [%s]: About to %s dataset %s in %s.\n",
	     __FILE__, xferTypeText, datasetName, fileTypeText);
    }
  }

  io_ncmpi_name_to_id_c(fileID, datasetName, &varID);

  /* Construct arrays for ncmpi_put_vars_all. */
  for (i=0; i<numDataDims; ++i) {
    if (globalOffset[i] > 0 && localSubSize[i] == 0) {
      /* Setting startDisk[i] to 0 stops an overly cautious index out of
	 bounds assertion error in pnetcdf library >= version 1.2.0pre1. */
      startDisk[i] = 0;
    } else {
      /* Global offset on disk */
      startDisk[i] = (MPI_Offset) globalOffset[i];
    }
    /* Amount of data myPE writes to disk */
    count[i] = (MPI_Offset) localSubSize[i];
  }

  /* We use a global MPI datatype (pointed to by pMPI_DataType) to
     pick up block data at location pData.  Make sure datatype is appropriate!
     localSubSize[0] contains the number of blocks on a given processor.
     startDisk, count refers to data in the file.
     pData, localSubSize[0], *pMPI_DataType refers to data in memory. */
  pMPI_DataType = io_get_grid_mpi_type(fileFmt, fileType, gridDataStruct);
  assert(pMPI_DataType != NULL);

  if (nonBlocking) {
    pNCMPI_Request = io_ncmpi_nonblocking_get_request();
    assert(pNCMPI_Request != NULL);

    /* The functions are macros that select the appropriate ncmpi
       iput/iget function */
    if (xferType == IO_WRITE_XFER) {
      err = flash_iput_vara(fileID, varID, startDisk, count, pData,
			    localSubSize[0], *pMPI_DataType,
			    pNCMPI_Request);
    } else if (xferType == IO_READ_XFER) {
      err = flash_iget_vara(fileID, varID, startDisk, count, pData,
			    localSubSize[0], *pMPI_DataType,
			    pNCMPI_Request);
    }
  } else {
    if (xferType == IO_WRITE_XFER) {
      err = ncmpi_put_vara_all(fileID, varID, startDisk, count, pData,
			       localSubSize[0], *pMPI_DataType);
    } else if (xferType == IO_READ_XFER) {
      err = ncmpi_get_vara_all(fileID, varID, startDisk, count, pData,
			       localSubSize[0], *pMPI_DataType);
    }
  }
  assert(err == NC_NOERR);
}
