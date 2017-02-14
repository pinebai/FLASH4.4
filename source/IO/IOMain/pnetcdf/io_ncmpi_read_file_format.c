#include "io_ncmpi_read_file_format.h"
#define DEBUG_IO

void FTOC(io_ncmpi_read_file_format)(const int * const pMyPE,
				     const int * const pFileID,
				     int *pFileFormat)
{
  const int myPE = *pMyPE;
  const int fileID = *pFileID;
  int err;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  err = ncmpi_get_att_int(fileID, NC_GLOBAL, "file_format_version", pFileFormat);
  assert (err == NC_NOERR);

  if (debugIO) {
    if (myPE == MASTER_PE) {
      printf(" [%s]: File format version is %d.\n",
	     __FILE__, *pFileFormat);
    }
  }
}
