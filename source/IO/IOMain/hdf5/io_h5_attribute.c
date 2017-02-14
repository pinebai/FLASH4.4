#include "io_h5_attribute.h"

/* WARNING: The caller must null-terminate all strings */

void io_h5_attribute_create(const int myPE,
			    const int fileID,
			    const int diskType,
			    const int dims,
			    const int diskSize[],
			    const char datasetName[],
			    const char attDatasetName[])
{
  const hid_t hFileID = (hid_t) fileID;
  hsize_t hDiskSize[IO_MAX_DIMS];
  hid_t hDiskType, dsetID, attID, attDataspace;
  herr_t err;
  int trueDims, parallelIO, i;

#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif

  if (parallelIO || (!parallelIO && myPE == MASTER_PE)) {
    assert(dims > 0 && dims <= IO_MAX_DIMS);
    for (i=0; i<dims; ++i) {
      hDiskSize[i] = (hsize_t) diskSize[i];
    }

    trueDims = dims;
    if (diskType == IO_FLASH_STRING) {
      hDiskType = io_h5_type_create_string(diskSize[dims-1]);
      if (dims == 1) {
	/* There is only 1 string */
	hDiskSize[0] = 1;
      } else if (dims > 1) {
	trueDims = trueDims - 1;
      }
    } else {
      /* We are using simple primitive types */
      hDiskType = io_h5_type_hid_primitive(diskType);
    }

    /* Open appropriate dataset and then add the attribute */
    dsetID = H5Dopen(hFileID, datasetName);
    assert(dsetID >= 0);

    attDataspace = H5Screate_simple(trueDims, hDiskSize, NULL);
    assert(attDataspace >= 0);

    attID = H5Acreate(dsetID, attDatasetName, hDiskType,
		      attDataspace, H5P_DEFAULT);
    assert(attID >= 0);


    /* Free allocated space. */
    err = H5Sclose(attDataspace);
    assert(err >= 0);

    err = H5Aclose(attID);
    assert(err >= 0);

    err = H5Dclose(dsetID);
    assert(err >= 0);

    if (diskType == IO_FLASH_STRING) {
      io_h5_type_free_string(hDiskType);
    }
  }
}


/* We pass the memory datatype so that the HDF5 library can
   convert between double in memory and float in file */
void io_h5_attribute_write(const int myPE,
			   const int fileID,
			   const int memType,
			   const char datasetName[],
			   const char attDatasetName[],
			   const void * const pData)
{
  const hid_t hFileID = (hid_t) fileID;
  hid_t hMemType, dsetID, attID;
  herr_t err;
  int parallelIO;

#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif

  if (parallelIO || (!parallelIO && myPE == MASTER_PE)) {

    /* Open appropriate dataset and then attribute */
    dsetID = H5Dopen(hFileID, datasetName);
    assert(dsetID >= 0);


#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8)
    /* This function is now deprecated by HDF5 */
    attID = H5Aopen_name(dsetID, attDatasetName);
#else
    attID = H5Aopen(dsetID, attDatasetName, H5P_DEFAULT);
#endif
    assert(attID >= 0);


    if (memType == IO_FLASH_STRING) {
      hMemType = H5Aget_type(attID);
      assert(hMemType >= 0);
    } else {
      hMemType = io_h5_type_hid_primitive(memType);
    }

    /* Write the attribute */
    err = H5Awrite(attID, hMemType, pData);
    assert(err >= 0);


    if (memType == IO_FLASH_STRING) {
      io_h5_type_free_string(hMemType);
    }

    err = H5Aclose(attID);
    assert(err >= 0);

    err = H5Dclose(dsetID);
    assert(err >= 0);
  }
}
