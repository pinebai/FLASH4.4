#include <assert.h>
#include <string.h>
#include <stdio.h>

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

void FTOC(io_ncmpi_create_dimids)(const int * const pFileID,
				  const int * const pDimVal,
				  const char dimName[],
				  const int * const pDimNameLen)
{
  const int fileID = *pFileID;
  const int dimVal = *pDimVal;
  const int dimNameLen = *pDimNameLen;
  int tmpVar, err;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif
  MPI_Offset mpiVal;
  char nullTermDimName[MAX_STRING_LENGTH+1];

  mpiVal = dimVal; /* Caller must be aware of possible int overflows */

  assert(dimNameLen > 0 && dimNameLen <= MAX_STRING_LENGTH);
  strncpy(nullTermDimName, dimName, dimNameLen);
  nullTermDimName[dimNameLen] = '\0';

  if (debugIO) {
    printf(" [%s]: Creating dimension id: %s with value: %d.\n",
           __FILE__, nullTermDimName, dimVal);
  }
  err = ncmpi_def_dim(fileID, nullTermDimName, mpiVal, &tmpVar);
  assert(err == NC_NOERR);
}
