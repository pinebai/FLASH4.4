#include "io_attribute.h"

/*
   We require the calling code to pass the length of each string
   as extra arguments.  We allow lengths up to MAX_STRING_LENGTH
   characters (excluding the null).

   The function prototype is chosen to be compatible with traditional
   Fortran (without interfaces) and Fortran 2003 (with interfaces).

   We ignore the hidden string length arguments that are typically
   passed by value at the end of the argument list.

   Portability is only guaranteed if the Fortran 2003 interfaces are used!
*/

void FTOC(io_attribute_create)(const int * const pMyPE,
			       const int * const pFileID,
			       const int * const pLibType,
			       const int * const pDiskType,
			       const int * const pDims,
			       const int datasetSize[],
			       const char datasetName[],
			       const int * const pDsetNameLen,
			       const char attDatasetName[],
  			       const int * const pAttNameLen)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int libType = *pLibType;
  const int diskType = *pDiskType;
  const int dims = *pDims;
  const size_t dsetNameLen = (size_t)*pDsetNameLen;
  const size_t attNameLen = (size_t)*pAttNameLen;
  char nullTermDsetName[MAX_STRING_LENGTH+1],
    nullTermAttName[MAX_STRING_LENGTH+1];

  assert(dsetNameLen > 0 && dsetNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermDsetName, datasetName, dsetNameLen);
  nullTermDsetName[dsetNameLen] = '\0';

  assert(attNameLen > 0 && attNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermAttName, attDatasetName, attNameLen);
  nullTermAttName[attNameLen] = '\0';



  if (libType == IO_FILE_PNETCDF) {

#if defined(FLASH_IO_PNETCDF)
    io_ncmpi_attribute_create(myPE, fileID, diskType, dims,
			      datasetSize, nullTermDsetName,
			      nullTermAttName);
#endif

  } else if (libType == IO_FILE_HDF5) {

#if defined(FLASH_IO_HDF5)
    io_h5_attribute_create(myPE, fileID, diskType, dims,
			   datasetSize, nullTermDsetName,
			   nullTermAttName);
#endif

  } else {
    Driver_abortFlashC("[io_attribute_create]: Unknown I/O");
  }
}


void FTOC(io_attribute_write)(const int * const pMyPE,
			      const int * const pFileID,
			      const int * const pLibType,
			      const int * const pMemType,
			      const char datasetName[],
			      const int * const pDsetNameLen,
			      const char attDatasetName[],
			      const int * const pAttNameLen,
			      const void * const pData)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  const int libType = *pLibType;
  const int memType = *pMemType;
  const size_t dsetNameLen = (size_t)*pDsetNameLen;
  const size_t attNameLen = (size_t)*pAttNameLen;
  char nullTermDsetName[MAX_STRING_LENGTH+1],
    nullTermAttName[MAX_STRING_LENGTH+1];

  assert(dsetNameLen > 0 && dsetNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermDsetName, datasetName, dsetNameLen);
  nullTermDsetName[dsetNameLen] = '\0';

  assert(attNameLen > 0 && attNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermAttName, attDatasetName, attNameLen);
  nullTermAttName[attNameLen] = '\0';

  if (libType == IO_FILE_PNETCDF) {

#if defined(FLASH_IO_PNETCDF)
    io_ncmpi_attribute_write(myPE, fileID, memType,
			     nullTermDsetName, nullTermAttName, pData);
#endif

  } else if (libType == IO_FILE_HDF5) {

#if defined(FLASH_IO_HDF5)
    io_h5_attribute_write(myPE, fileID, memType,
			  nullTermDsetName, nullTermAttName, pData);
#endif

  } else {
    Driver_abortFlashC("[io_attribute_write]: Unknown I/O");
  }
}
