#include <stdio.h>
#include "mangle_names.h"
#include "constants.h"
#include <hdf5.h>


void FTOC(io_h5create_raydset)(const int * const pFileID)
{
  /* Create the data space with unlimited dimensions. */
  hsize_t dims[2] = {0, 5};
  hsize_t maxdims[2] = {H5S_UNLIMITED, 5};
  hid_t dataspace = H5Screate_simple (2, dims, maxdims);
  
  /* Modify dataset creation properties, i.e. enable chunking  */
  hsize_t chunk_dims[2] = {256, 5}; /* TODO: Pass in good chunk size estimate */
  hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk ( cparms, 2, chunk_dims);
  
  /* Create a new dataset within the file using cparms
     creation properties.  */
  hid_t dsetid = H5Dcreate (*pFileID, "RayData", H5T_NATIVE_DOUBLE, dataspace,
			    cparms);

  H5Sclose(dataspace);
  H5Dclose(dsetid);
}
