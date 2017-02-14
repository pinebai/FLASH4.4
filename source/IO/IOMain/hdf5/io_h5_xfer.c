#include "io_h5_xfer.h"

int io_h5_xfer(const int myPE, const hid_t fileID, const int xferType,
	       const hid_t hXferList, const char datasetName[],
	       const hid_t hMemType, const hsize_t hMemSize[],
	       const hsize_t hMemStart[], const hsize_t hMemCount[],
	       const hsize_t hDiskStart[], const hsize_t hDiskCount[],
	       const int dims, void * pData)
{
  herr_t (*old_func)(void*) = NULL;
  void *old_client_data = NULL;
  hsize_t hDummyMemSize[IO_MAX_DIMS];
  hid_t dataspace, memspace, dataset;
  herr_t err;
  int i, zeroSize = 0;
  int ierr;


  for (i=0; i<dims; ++i) {
    /* This function only deals with primitive types so following assert
       must be true */
    assert (hMemCount[i] == hDiskCount[i]);
    if (hMemCount[i] == 0 || hDiskCount[i] == 0) zeroSize = 1;
  }


  /* open the dataset
     suppress the error message if the variable does not exist in file */
  err = H5Eget_auto(&old_func, &old_client_data);
  assert(err >= 0);
  err = H5Eset_auto(NULL, NULL);
  assert(err >= 0);
  dataset = H5Dopen(fileID, datasetName);
  err = H5Eset_auto(old_func, old_client_data);
  assert(err >= 0);


  if (dataset < 0) {

    if (xferType == IO_READ_XFER) {
      ierr = -1;
      if (myPE == MASTER_PE) {
	printf(" [%s]: Skipping missing dataset '%s'.\n",
	       __FILE__, datasetName);
      }
    } else {
      printf("[%s]: Processor %d failed during H5Dopen on dataset %s.\n",
	     __FILE__, myPE, datasetName);
      ierr = Driver_abortFlashC("Error! H5Dopen failed");
    }

  } else {

    dataspace = H5Dget_space(dataset);
    assert (dataspace >= 0);


    for (i=0; i<dims; ++i) {
      if (zeroSize) {
	hDummyMemSize[i] = 1; /* Must be > 0 to satisfy H5Screate_simple */
      } else {
	hDummyMemSize[i] = hMemSize[i];
      }
    }
    memspace = H5Screate_simple(dims, hDummyMemSize, NULL);
    assert (memspace >= 0);


    if (zeroSize) {
      /* Select nothing if myPE is not contributing any data to disk. */
      err = H5Sselect_none(dataspace);
      assert (err >= 0);
      err = H5Sselect_none(memspace);
      assert (err >= 0);
    } else {
      /* Define myPE's portion in global data */
      err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hDiskStart,
				NULL, hDiskCount, NULL);
      assert (err >= 0);

      /* Define how myPE's data is stored in memory */
      err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
				hMemStart, NULL, hMemCount, NULL);
      assert (err >= 0);
    }


    if (xferType == IO_READ_XFER) {
      /* Read data from disk hyperslab to memory */
      err = H5Dread(dataset, hMemType, memspace, dataspace,
		    hXferList, pData);
    } else if (xferType == IO_WRITE_XFER) {
      /* Write data from memory to disk hyperslab */
      err = H5Dwrite(dataset, hMemType, memspace, dataspace,
		     hXferList, pData);
    }
    assert (err >= 0);


    err = H5Sclose(dataspace);
    assert (err >= 0);

    err = H5Dclose(dataset);
    assert (err >= 0);

    err = H5Sclose(memspace);
    assert (err >= 0);

    ierr = 0;
  }

  return ierr;
}
