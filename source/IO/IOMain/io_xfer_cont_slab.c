#include "io_xfer_cont_slab.h"

/* Write a fully-contiguous section */
void FTOC(io_xfer_cont_slab)(const int * const pMyPE,
			     const int * const pFileID,
			     const int * const pLibType,
			     const int * const pXferType,
			     const int * const pTypeMatchedXfer,
			     const char datasetName[],
			     const int * const pNameLength,
			     const int * const pMemType,
			     const int memSize[],
			     const int memStart[],
			     const int memCount[],
			     const int diskStart[],
			     const int diskCount[],
			     const int * const pDims,
			     void * pData,
			     int * const pErr)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int libType = *pLibType;
  const size_t nameLength = (size_t)*pNameLength;
  const int memType = *pMemType;
  const int dims = *pDims;
  const int xferType = *pXferType;
  const int typeMatchedXfer = *pTypeMatchedXfer;
  int i, err;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif
  char nullTermDatasetName[MAX_STRING_LENGTH+1];


  assert (nameLength > 0 && nameLength <= MAX_STRING_LENGTH);
  strncpy (nullTermDatasetName, datasetName, nameLength);
  nullTermDatasetName[nameLength] = '\0';


  for (i=0; i<dims; ++i) {
    if (debugIO) {
      printf(" [io_xfer_cont_slab]: Dataset %s. Proc %d. Dim %d. "  \
               "Disk (start %d, count %d). "                   \
	     "Mem (size %d, start %d, count %d).\n",
	     nullTermDatasetName, myPE, i,
	     diskStart[i], diskCount[i],
	     memSize[i], memStart[i], memCount[i]);
    }

    assert (diskStart[i] >= 0);
    assert (diskCount[i] >= 0);
    assert (memSize[i] >= 0);
    assert (memStart[i] >= 0);
    assert (memCount[i] >= 0);
    assert (memSize[i] >= memStart[i] + memCount[i]);
  }

  if (libType == IO_FILE_PNETCDF) {

#if defined (FLASH_IO_PNETCDF)
    err = io_ncmpi_xfer_wrapper(myPE, fileID, xferType, nullTermDatasetName,
				memType, memSize, memStart, memCount,
				diskStart, diskCount, dims, pData);
#endif

  } else if (libType == IO_FILE_HDF5) {

#if defined(FLASH_IO_HDF5)
    err = io_h5_xfer_wrapper(myPE, fileID, xferType, typeMatchedXfer,
			     nullTermDatasetName,
			     memType, memSize, memStart, memCount,
			     diskStart, diskCount, dims, pData);
#endif

  } else {
    Driver_abortFlashC("[io_xfer_cont_slab]: Unknown I/O");
  }
  *pErr = err;
}
