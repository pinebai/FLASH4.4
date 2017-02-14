#include <hdf5.h>
#include <assert.h>
#include "mangle_names.h"
#include "constants.h"
#include "Flash.h"

void FTOC(io_h5read_present_dims)(const int * const pFileID,
				  int *pDims)
{
  hsize_t dimens_2d[2], maxdimens_2d[2];
  hid_t fileID, dataspace, dataset;
  herr_t herr;
  int ierr;

  assert (sizeof(int) == sizeof(hid_t));
  fileID = (hid_t) *pFileID;

  dataset = H5Dopen(fileID, "coordinates");
  assert (dataset >= 0);

  dataspace = H5Dget_space(dataset);
  assert (dataspace >= 0);

  ierr = H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
  assert (ierr >= 0);

  herr = H5Sclose(dataspace);
  assert (herr >= 0);

  herr = H5Dclose(dataset);
  assert (herr >= 0);

  *pDims = dimens_2d[1];
}
