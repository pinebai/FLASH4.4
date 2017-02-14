#include "io_h5_read_file_format.h"
#define DEBUG_IO

void FTOC(io_h5_read_file_format)(const int * const pMyPE,
				  const int * const pFileID,
				  int *pFileFormat)
{
  const int myPE = *pMyPE;
  const hid_t fileID = (hid_t) *pFileID;
  hid_t dset, dtypeInt;
  int err;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  dset = H5Dopen(fileID, "sim info");
  assert(dset >= 0);
  dtypeInt = H5Tcreate(H5T_COMPOUND, sizeof(int));
  assert(dtypeInt >= 0);
  err = H5Tinsert(dtypeInt, "file format version", 0, H5T_NATIVE_INT);
  assert(err >= 0);


  err = H5Dread(dset, dtypeInt, H5S_ALL, H5S_ALL, H5P_DEFAULT, pFileFormat);
  assert(err >= 0);


  err = H5Tclose(dtypeInt);
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);

  if (debugIO) {
    if (myPE == MASTER_PE) {
      printf(" [%s]: File format version is %d.\n",
	     __FILE__, *pFileFormat);
    }
  }
}
