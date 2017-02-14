#include "io_ncmpi_define_mode.h"

void FTOC(io_ncmpi_define_mode_enddef)(const int * const pFileID)
{
  const int fileID = *pFileID;
  int err;

  err = ncmpi_enddef(fileID);
  assert(err == NC_NOERR);
}

void FTOC(io_ncmpi_define_mode_redef)(const int * const pFileID)
{
  const int fileID = *pFileID;
  int err;

  err = ncmpi_redef(fileID);
  assert(err == NC_NOERR);
}
