#include "io_h5_type_matched_xfer.h"

/* This function works around limitations of the HDF5 library.  We have
   noticed that HDF5 resorts to independent parallel I/O when the memory
   type is different from the disk type.  We therefore allocate and use a
   temporary buffer of the same type as the disk type so that HDF5
   can perform the data transfer using collective I/O optimizations. */

int
io_h5_type_matched_xfer(const int myPE, const hid_t hFileID, const int xferType,
			const hid_t hXferList, const char datasetName[],
			const hid_t hMemType, const hsize_t hMemSize[],
			const hsize_t hMemStart[], const hsize_t hMemCount[],
			const hsize_t hDiskStart[], const hsize_t hDiskCount[],
			const int dims, void * pData)
{
  hid_t hMatchedMemType;
  hsize_t hPackedSize;
  herr_t herr;
  void *pHdf5Buf = NULL;
  double *pTmpDblBuf = NULL;
  float *pTmpFltBuf = NULL;
  int useExtraFloatBuffer, err, packedSize, i;

  MPI_Datatype mpiMemTypeScalar = MPI_DOUBLE, mpiMemType;
  int memSize[IO_MAX_DIMS], memCount[IO_MAX_DIMS], memStart[IO_MAX_DIMS];
  int order = MPI_ORDER_C, one = 1;

  assert(dims <= IO_MAX_DIMS);
  for (i=0; i<dims; ++i) {
    assert(hMemSize[i] >= 0);
    assert(hMemStart[i] >= 0);
    assert(hMemCount[i] >= 0);
    assert(hMemSize[i] >= hMemCount[i]);
    assert(hMemCount[i] == hDiskCount[i]);
  }


  useExtraFloatBuffer = io_h5_use_extra_float_buffer(myPE, hFileID,
						     datasetName, hMemType);

  /* Pre-transfer */
  if (useExtraFloatBuffer) {

    /* Create a MPI subarray to describe the requested HDF5 hyperslab */
    hPackedSize = 1;
    for (i=0; i<dims; ++i) {
      hPackedSize *= hMemCount[i];
    }
    assert(hPackedSize <= INT_MAX);
    packedSize = (int) hPackedSize;

    if (packedSize >= 1) {
      for (i=0; i<dims; ++i) {
	assert(hMemSize[i] <= INT_MAX);
	assert(hMemStart[i] <= INT_MAX);
	assert(hMemCount[i] <= INT_MAX);
	memSize[i] = (int) hMemSize[i];
	memStart[i] = (int) hMemStart[i];
	memCount[i] = (int) hMemCount[i];
      }
      err = MPI_Type_create_subarray(dims, memSize, memCount,
				     memStart, order, mpiMemTypeScalar,
				     &mpiMemType);
      assert(err == MPI_SUCCESS);
      err = MPI_Type_commit(&mpiMemType);
      assert(err == MPI_SUCCESS);

      pTmpDblBuf = malloc (packedSize * sizeof(*pTmpDblBuf));
      assert(NULL != pTmpDblBuf);
      pTmpFltBuf = malloc (packedSize * sizeof(*pTmpFltBuf));
      assert(NULL != pTmpFltBuf);
      if (xferType == IO_WRITE_XFER) {
	/* Pack the data into a contiguous buffer of type float */
	err = io_repack_data(pData, one, mpiMemType, (void*) pTmpDblBuf,
			     packedSize, mpiMemTypeScalar);
	assert(err == 0);
	for (i=0; i<packedSize; ++i) {
	  pTmpFltBuf[i] = (float) pTmpDblBuf[i];
	}
      }
    } else {
      /* There is no need to commit basic datatypes -
	 they are "pre-committed". */
      mpiMemType = MPI_DATATYPE_NULL;
    }

    hMatchedMemType = H5Tcopy(H5T_NATIVE_FLOAT);
    assert(hMatchedMemType >= 0);
    pHdf5Buf = pTmpFltBuf;

  } else {

    hMatchedMemType = H5Tcopy(hMemType);
    assert(hMatchedMemType >= 0);
    pHdf5Buf = pData;
  }


  /* Transfer */
  err = io_h5_xfer(myPE, hFileID, xferType, hXferList, datasetName,
		   hMatchedMemType, hMemSize, hMemStart, hMemCount,
		   hDiskStart, hDiskCount, dims, pHdf5Buf);
  assert(err == 0);


  /* Post-transfer */
  if (useExtraFloatBuffer) {
    if (packedSize >= 1) {

      if (xferType == IO_READ_XFER) {
	/* Unpack the data from a contiguous buffer of type float */
	for (i=0; i<packedSize; ++i) {
	  pTmpDblBuf[i] = (double) pTmpFltBuf[i];
	}
	err = io_repack_data((void*) pTmpDblBuf, packedSize, mpiMemTypeScalar,
			     pData, one, mpiMemType);
	assert(err == 0);
      }

      err = MPI_Type_free(&mpiMemType);
      assert(err == MPI_SUCCESS);
      free(pTmpDblBuf);
      free(pTmpFltBuf);
    }
  }


  herr = H5Tclose(hMatchedMemType);
  assert(herr >= 0);

  pTmpFltBuf = NULL;
  pTmpDblBuf = NULL;
  pHdf5Buf = NULL;

  return 0;
}



/* The purpose of this code is to find out whether the memory type
   is double and the disk type is float. */

int
io_h5_use_extra_float_buffer(const int myPE, const hid_t hFileID,
			     const char datasetName[], const hid_t hMemType)
{
  htri_t cmp1, cmp2, cmp3;
  hid_t dataset, hDiskType, hDiskTypeNative;
  int useExtraFloatBuffer;
  herr_t herr;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  useExtraFloatBuffer = 0;


  dataset = H5Dopen(hFileID, datasetName);
  assert (dataset >= 0);
  hDiskType = H5Dget_type(dataset);
  assert(hDiskType >= 0);
  hDiskTypeNative = H5Tget_native_type(hDiskType, H5T_DIR_ASCEND);
  assert(hDiskTypeNative >= 0);

  cmp1 = H5Tequal(hMemType,hDiskTypeNative);
  assert(cmp1 >= 0);
  if (cmp1 == 0) {
    if (debugIO && myPE == MASTER_PE) {
      printf(" [%s]: Dataset %s (Memory type differs from file type)\n",
	     __FILE__, datasetName);
    }
    cmp2 = H5Tequal(hMemType,H5T_NATIVE_DOUBLE);
    assert(cmp2 >= 0);
    cmp3 = H5Tequal(hDiskTypeNative,H5T_NATIVE_FLOAT);
    assert(cmp3 >= 0);
    if (cmp2 > 0 && cmp3 > 0) {
      useExtraFloatBuffer = 1;
      if (debugIO && myPE == MASTER_PE) {
	printf(" [%s]: Dataset %s (Using extra float buffer to enable "
	       "collective I/O)\n", __FILE__, datasetName);
      }
    }
  }

  herr = H5Tclose(hDiskTypeNative);
  assert(herr >= 0);
  herr = H5Tclose(hDiskType);
  assert(herr >= 0);
  herr = H5Dclose(dataset);
  assert(herr >= 0);


  return useExtraFloatBuffer;
}
