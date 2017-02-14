#include "io_ncmpi_xfer.h"

int io_ncmpi_xfer(const int myPE, const int fileID, const int xferType,
		  const int varID, const MPI_Datatype memType,
		  const MPI_Offset memCountScalar, const MPI_Offset diskStart[],
		  const MPI_Offset diskCount[], const int dims, void * pData)
{
  int err;

  if (xferType == IO_READ_XFER) {
    /* Read data from file to memory */
    err = ncmpi_get_vara_all(fileID, varID, diskStart, diskCount,
			     pData, memCountScalar, memType);
  } else if (xferType == IO_WRITE_XFER) {
    /* Write data from memory to file */
    err = ncmpi_put_vara_all(fileID, varID, diskStart, diskCount,
			     pData, memCountScalar, memType);
  }
  if (err != NC_NOERR) printf("%s\n", ncmpi_strerror(err));
  assert (err == NC_NOERR);
  return 0;
}
