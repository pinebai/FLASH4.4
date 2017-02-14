#include "io_h5_xfer_wrapper.h"

int io_h5_xfer_wrapper(const int myPE, const int fileID, const int xferType,
		       const int typeMatchedXfer,
		       const char datasetName[], const int memType,
		       const int memSize[], const int memStart[],
		       const int memCount[], const int diskStart[],
		       const int diskCount[], const int dims, void * pData)
{
  hsize_t hDiskStart[IO_MAX_DIMS], hDiskCount[IO_MAX_DIMS],
    hMemSize[IO_MAX_DIMS], hMemStart[IO_MAX_DIMS], hMemCount[IO_MAX_DIMS],
    hStrLen;
  hid_t hXferList, hMemType, hFileID;
  herr_t hErr;
  int i, j, err, errRtn, totStr, parallelIO, xferTypeLowLevel, rank;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif
  char *pChar = NULL;
  char nullTermDebugStr[1000];


  /* We make hid_t == int assumption everywhere in FLASH.  This is
     safe because hid_t is a typedef of int in HDF5.  The assert
     will notify us of possible future changes to hid_t */
  assert (sizeof(int) == sizeof(hid_t));
  hFileID = (hid_t) fileID;


  /* Check if only MASTER_PE used for read/write.  If so ensure we do not
     synchronise across all processors, i.e. set parallelIO = 0. */
  if (xferType == IO_READ_XFER_MASTER_PE) {
    xferTypeLowLevel = IO_READ_XFER;
    parallelIO = 0;
  } else if (xferType == IO_WRITE_XFER_MASTER_PE) {
    xferTypeLowLevel = IO_WRITE_XFER;
    parallelIO = 0;
  } else {
    xferTypeLowLevel = xferType;
#if defined(IO_HDF5_PARALLEL)
    parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
    parallelIO = 0;
#else
    parallelIO = 0;
#endif
  }


  /* Return from this function if we are in serialIO mode AND we are
     not the master processor */
  errRtn = 0;
  if (parallelIO || (!parallelIO && myPE == MASTER_PE)) {
    assert (myPE >= 0);
    assert (dims >= 1 && dims <= IO_MAX_DIMS);
    rank = dims;

    /* Create arrays of type hsize_t describing myPE's portion in global data */
    for (i=0; i<rank; ++i) {
      hDiskStart[i] = (hsize_t) diskStart[i];
      hDiskCount[i] = (hsize_t) diskCount[i];
      hMemSize[i] = (hsize_t) memSize[i];
      hMemStart[i] = (hsize_t) memStart[i];
      hMemCount[i] = (hsize_t) memCount[i];
    }

    /* Set independent/collective transfer mode */
    hXferList = H5Pcreate(H5P_DATASET_XFER);
    assert (hXferList != -1);
#ifdef H5_HAVE_PARALLEL
    if (parallelIO && HDF5_MODE == COLLECTIVE) {
      hErr = H5Pset_dxpl_mpio(hXferList, H5FD_MPIO_COLLECTIVE);
      assert (hErr >= 0);
    } else {
      hErr = H5Pset_dxpl_mpio(hXferList, H5FD_MPIO_INDEPENDENT);
      assert (hErr >= 0);
    }
#endif

    if (memType == IO_FLASH_STRING) {
      for (i=0; i<rank-1; ++i) {
	/* Assumptions that are only needed for string types */
	assert (memSize[i] == memCount[i]);
	assert (diskCount[i] == memCount[i]);
      }
      if (debugIO) {
	totStr = 1;
	pChar = pData;
	for (i=0; i<rank-1; ++i) {
	  totStr *= memCount[i];
	}
	for (j=0; j<totStr; ++j) {
	  /* Contiguous dimension is the length of each string */
	  strncpy(nullTermDebugStr, &pChar[(j*memCount[rank-1])],
		  (size_t)memCount[rank-1]);
	  nullTermDebugStr[memCount[rank-1]] = '\0';
	  if (myPE == MASTER_PE) {
	    printf(" [%s]: %s (string %d) is %s\n",
		   __FILE__, datasetName, j, nullTermDebugStr);
	  }
	}
      }

      hMemType = io_h5_type_create_string(memCount[rank-1]);
      hStrLen = hMemCount[rank-1];
      if (rank == 1) {
	/* There is only 1 string of length hStrLen */
	hMemSize[0] = 1; hMemStart[0] = 0;
	hMemCount[0] = 1; hDiskCount[0] = 1;
	hDiskStart[0] = hDiskStart[0] / hStrLen;
      } else if (rank > 1) {
        rank = rank - 1;
      }
    } else {
      hMemType = io_h5_type_hid_primitive(memType);
    }


    /* This is changeable at runtime.  It is set to 1 (i.e. true) by default
       so that we get the full benefits of collective I/O at scale despite the
       limitations of HDF5 library. */
    if (typeMatchedXfer) {
      err = io_h5_type_matched_xfer(myPE, hFileID, xferTypeLowLevel, hXferList,
				    datasetName,
				    hMemType, hMemSize, hMemStart, hMemCount,
				    hDiskStart, hDiskCount, rank, pData);
      assert (err == 0);
    }
    else {
      errRtn = io_h5_xfer(myPE, hFileID, xferTypeLowLevel, hXferList, datasetName,
		       hMemType, hMemSize, hMemStart, hMemCount,
		       hDiskStart, hDiskCount, rank, pData);
    }
    if (debugIO) {
      err = io_h5_report_xfer_method(myPE, hXferList, datasetName);
      assert (err == 0);
    }

    hErr = H5Pclose(hXferList);
    assert (hErr >= 0);
    if (memType == IO_FLASH_STRING) {
      io_h5_type_free_string(hMemType);
    }
  } /* if (parallelIO || (!parallelIO && myPE == MASTER_PE)) */

  return errRtn;
}
