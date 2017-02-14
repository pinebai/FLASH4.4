#include "io_xfer_mesh_dataset.h"

/* pData is void * to be compatible with F2003 Fortran / C 
interoperability requirements.  In Fortran we are specifically 
passing a memory address rather than an array. */

void FTOC(io_xfer_mesh_dataset)(const int * const pMyPE,
				const int * const pFileID,
				const int * const pLibType,
				const int * const pXferType,
				const int * const pFileType,
				const int * const pFileFmt,
				const int * const pGridDataStruct,
				const int * const pNumDataDims,
				const int * const pNumGridVars,
				const int * const pNonBlocking,
				const int * const pPrePackData,
				const int diskOffset[],
				const int memSize[],
				const int memSubSize[],
				const int memOffset[],
				const int memVarOffset[],
				const char datasetName[],
				const int * const pDsetNameLen,
				void * pData)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int libType = *pLibType;
  const int xferType = *pXferType;
  const int fileType = *pFileType;
  const int fileFmt = *pFileFmt;
  const int gridDataStruct = *pGridDataStruct;
  const int numDataDims = *pNumDataDims;
  const int numGridVars = *pNumGridVars;
  const int nonBlocking = *pNonBlocking;
  const int prePackData = *pPrePackData;
  const int dsetNameLen = *pDsetNameLen;
  double *pDblData = pData;
  char nullTermDsetName[MAX_STRING_LENGTH+1];

  assert(dsetNameLen > 0 && dsetNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermDsetName, datasetName, dsetNameLen);
  nullTermDsetName[dsetNameLen] = '\0';
  

  if (libType == IO_FILE_PNETCDF) {

#if defined(FLASH_IO_PNETCDF)
    io_ncmpi_xfer_mesh_dataset(myPE, fileID, xferType, fileType,
			       fileFmt, gridDataStruct, numDataDims,
			       nonBlocking, diskOffset, memSubSize,
			       nullTermDsetName, pDblData);
#endif


  } else if (libType == IO_FILE_HDF5) {
    if (prePackData == 1) {

#if defined(FLASH_IO_HDF5)
      io_h5_xfer_packed_mesh_dataset(myPE, fileID, xferType, fileFmt,
				     fileType, gridDataStruct,
				     numDataDims, nullTermDsetName,
				     diskOffset, memSubSize, pDblData);
#endif

    } else {

#if defined(FLASH_IO_HDF5)
      io_h5_xfer_mesh_dataset(myPE, fileID, xferType, fileFmt,
			       fileType, gridDataStruct, numDataDims,
			       numGridVars, memVarOffset,
			       diskOffset, memSize, memSubSize, memOffset,
			       nullTermDsetName, pDblData);
#endif
    } 
  } else {
    Driver_abortFlashC("[io_xfer_mesh_dataset]: Unknown I/O");
  }
}
