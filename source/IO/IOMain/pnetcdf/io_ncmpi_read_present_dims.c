#include <pnetcdf.h>
#include "mangle_names.h"
#include "constants.h"
#include "Flash.h"

void FTOC(io_ncmpi_read_present_dims)(const int * const pFileID,
				      int * pDims)
{
  const int fileID = *pFileID;
  int err, dims, dimID;

  /* Some of our old files used dim_NDIM as opposed to the current dim_MDIM.
     Attempting to read in too much data will cause an abort.*/
  err = ncmpi_inq_dimid(fileID, "dim_MDIM", &dimID);
  if (err == NC_NOERR) {
    dims = MDIM;
  } else {
    dims = NDIM;
  }
  *pDims = dims;
}
