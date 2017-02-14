#include "io_create_dataset.h"

void FTOC(io_create_dataset)(const int * const pMyPE,
			     const int * const pFileID,
			     const int * const pLibType,
			     const int * const pDiskType,
			     const int * const pDims,
			     const int dimIDs[],
			     const char datasetName[],
			     const int * const pDsetNameLen)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int libType = *pLibType;
  const int diskType = *pDiskType;
  const int dims = *pDims;
  const size_t dsetNameLen = (size_t)*pDsetNameLen;
  char nullTermDsetName[MAX_STRING_LENGTH+1];

  assert(dsetNameLen > 0 && dsetNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermDsetName, datasetName, dsetNameLen);
  nullTermDsetName[dsetNameLen] = '\0';


  if (libType == IO_FILE_PNETCDF) {

#if defined(FLASH_IO_PNETCDF)
    io_ncmpi_create_dataset(myPE, fileID, diskType, dims,
			    dimIDs, nullTermDsetName);
#endif

  } else if (libType == IO_FILE_HDF5) {

#if defined(FLASH_IO_HDF5)
    io_h5create_dataset(myPE, fileID, diskType, dims,
			dimIDs, nullTermDsetName);
#endif

  } else {
    Driver_abortFlashC("[io_create_dataset]: Unknown I/O");
  }
}
