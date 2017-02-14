#include "io_h5create_dataset.h"

void io_h5create_dataset(const int myPE,
			 const int fileID,
			 const int diskType,
			 const int dims,
			 const int diskSize[],
			 const char datasetName[])
{
  const hid_t hFileID = (hid_t) fileID;
  hsize_t hDiskSize[IO_MAX_DIMS];
  hid_t hDiskType, dataspace, dataset;
  herr_t err;
  int rank, i;

#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif


  assert (dims > 0);
  for (i=0; i<dims; ++i) {
    if (debugIO && myPE == MASTER_PE) {
      printf(" [io_h5create_dataset]: Dataset %s. Proc %d. Dim %d. Size %d.\n",
	     datasetName, myPE, i, diskSize[i]);
    }
    assert(diskSize[i] > 0);
    hDiskSize[i] = (hsize_t) diskSize[i];
  }
  rank = dims;


  if (diskType == IO_FLASH_STRING) {
    /* Create a HDF5 string type.  It is possible to have
       an N-dimensional array of strings (the final dimension
       contains the string length) */
    hDiskType = io_h5_type_create_string(diskSize[rank-1]);
    if (rank == 1) {
      /* There is only 1 string */
      hDiskSize[0] = 1;
    } else if (rank > 1) {
      rank = rank - 1;
    }
  } else {
    /* We are using simple primitive types */
    hDiskType = io_h5_type_hid_primitive(diskType);
  }


  /* Create the actual dataset (size & type in the file). */
  dataspace = H5Screate_simple(rank, hDiskSize, NULL);
  assert(dataspace >= 0);
  dataset = H5Dcreate(hFileID, datasetName, hDiskType,
                      dataspace, H5P_DEFAULT);
  assert(dataset >= 0);


  err = H5Sclose(dataspace);
  assert(err >= 0);
  err = H5Dclose(dataset);
  assert(err >= 0);

  if (diskType == IO_FLASH_STRING) {
    io_h5_type_free_string(hDiskType);
  }
}
