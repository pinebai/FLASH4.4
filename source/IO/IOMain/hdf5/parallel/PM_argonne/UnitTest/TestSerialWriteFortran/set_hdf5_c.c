#include "set_hdf5_c.h"

int HDF5_MODE = 0;

void FTOC(init_hdf5_c)(int *pFileID)
{
  const char filename[] = "Test_write_fortran.hdf5";
  *pFileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  assert(*pFileID >= 0);
}

void FTOC(finalise_hdf5_c)(int *pFileID)
{
  herr_t err;
  err = H5Fclose(*pFileID);
  assert(err >= 0);
}
