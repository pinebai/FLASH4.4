#include "io_h5_xfer_packed_mesh_dataset.h"

  /* This subroutine manually packs the data into a contiguous buffer
     so that it is a simple write operation for HDF5.  It requires
     allocating 2 or 3 bufffers.  If we are packing data into a contiguous
     buffer we cast it to the type that is stored in the file -- this
     allows collective HDF5 optimisations.  Whenever the memory type
     differs from the file type HDF5 falls back to independent mode

     We assume that data is stored as double in memory, but allow
     transfers to float or double in file.
*/

void io_h5_xfer_packed_mesh_dataset(const int myPE,
				    const int fileID,
				    const int xferType,
				    const int fileFmt,
				    const int fileType,
				    const int gridDataStruct,
				    const int numDataDims,
				    const char datasetName[],
				    const int globalOffset[],
				    const int localSubSize[],
				    double * pData)
{
  const hid_t hFileID = (hid_t) fileID;
  MPI_Datatype *pMPI_DataType = NULL;
  void *pHdf5Buf = NULL;
  /* pTmpDbleBuf and pTmpFltBuf must be set to NULL so that free() can be
     used even if the buffers are not allocated. */
  double *pTmpDbleBuf = NULL;
  float *pTmpFltBuf = NULL;
  herr_t err;
  size_t is, tmpBufSize;
  hid_t dataset, dtype, dtypeNative;
  int i, flashDiskType, typeMatchedXfer;
  int memOffset[IO_MESH_DIMS] = {0,0,0,0,0};
  char *fileTypeText = NULL, *xferTypeText = NULL;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  /* Set to 0 because this subroutine already handles the type matching */
  typeMatchedXfer = 0;

  /* Check input variables are sensible */
  assert(numDataDims <= IO_MESH_DIMS); /* Too many dims */
  assert(myPE >= 0); /* Invalid process identifier */
  assert(fileType == CHECKPOINTFILE || fileType == PLOTFILE);
  assert(gridDataStruct == CENTER || gridDataStruct == FACEX ||
	 gridDataStruct == FACEY || gridDataStruct == FACEZ ||
	 gridDataStruct == SCRATCH);

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


  /* Open the dataset so that we can obtain the datatype */
  dataset = H5Dopen(hFileID, datasetName);
  assert(dataset >= 0);
  dtype = H5Dget_type(dataset);
  assert(dtype >= 0);
  dtypeNative = H5Tget_native_type(dtype, H5T_DIR_ASCEND);
  assert(dtypeNative >= 0);

  flashDiskType = io_h5_type_flash_primitive(dtypeNative);

  err = H5Tclose(dtypeNative);
  assert(err >= 0);
  err = H5Tclose(dtype);
  assert(err >= 0);
  err = H5Dclose(dataset);
  assert(err >= 0);


  /* Select the datatype that describes data in memory to be
   read from / written to file */
  pMPI_DataType = io_get_grid_mpi_type(fileFmt, fileType, gridDataStruct);
  assert(NULL != pMPI_DataType);
  assert(MPI_DATATYPE_NULL != *pMPI_DataType);


  tmpBufSize = 1;
  for (i=0; i<numDataDims; ++i) {
    assert(localSubSize[i] >= 0);
    tmpBufSize *= localSubSize[i];
  }
  if (tmpBufSize > 0) {
    pTmpDbleBuf = malloc (tmpBufSize * sizeof(*pTmpDbleBuf));
    assert(NULL != pTmpDbleBuf);
  }


  /* The code here ensures that the primitive memory type matches the
     primitive file type */
  if (flashDiskType == IO_FLASH_DOUBLE) {
    if (debugIO && myPE == MASTER_PE) {
      printf(" [%s]: Dataset %s is double.\n", __FILE__, datasetName);
    }

    /* Double in memory and file so just point HDF5 buffer to the
       temporary buffer */
    pHdf5Buf = pTmpDbleBuf;

  } else if (flashDiskType == IO_FLASH_FLOAT) {
    if (debugIO && myPE == MASTER_PE) {
      printf(" [%s]: Dataset %s is float.\n", __FILE__, datasetName);
    }

    /* Need an extra buffer to perform a Double -> Float conversion */
    if (tmpBufSize > 0) {
      pTmpFltBuf = malloc (tmpBufSize * sizeof(*pTmpFltBuf));
      assert(NULL != pTmpFltBuf);
    }
    pHdf5Buf = pTmpFltBuf;
  }
  else {
    Driver_abortFlashC("File datatype is neither float nor double\n");
  }


  if (tmpBufSize > 0) {
    if (xferType == IO_WRITE_XFER) {
      /* Repack grid data into contiguous double precision data buffer */
      err = io_repack_data(pData, localSubSize[0], *pMPI_DataType,
			   pTmpDbleBuf, tmpBufSize, MPI_DOUBLE);
      if (flashDiskType == IO_FLASH_FLOAT) {
	for (is=0; is<tmpBufSize; ++is) {
	  pTmpFltBuf[is] = (float) pTmpDbleBuf[is];
	}
      }
    }
  }


  /* pHdf5Buf points to pTmpDbleBuf or pTmpFltBuf
     All MPI processes must call io_h5_xfer_wrapper */
  io_h5_xfer_wrapper(myPE, hFileID, xferType,
		     typeMatchedXfer,
		     datasetName, flashDiskType,
		     localSubSize, memOffset, localSubSize,
		     globalOffset, localSubSize,
		     numDataDims, pHdf5Buf);


  if (tmpBufSize > 0) {
    if (xferType == IO_READ_XFER) {
      if (flashDiskType == IO_FLASH_FLOAT) {
	for (is=0; is<tmpBufSize; ++is) {
	  pTmpDbleBuf[is] = (double) pTmpFltBuf[is];
	}
      }
      /* Repack grid data into mesh data structure */
      err = io_repack_data(pTmpDbleBuf, tmpBufSize, MPI_DOUBLE,
			   pData, localSubSize[0], *pMPI_DataType);
    }
  }


  free (pTmpDbleBuf);
  free (pTmpFltBuf);
  pHdf5Buf = NULL;
}
