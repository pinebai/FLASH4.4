#include "io_ncmpi_create_dataset.h"

void io_ncmpi_create_dataset(const int myPE,
			     const int fileID,
			     const int diskType,
			     const int dims,
			     const int dimIDs[],
			     const char datasetName[])
{
  int tmpVar, err;
  nc_type ncDataType;

  ncDataType = io_ncmpi_type_primitive(diskType);

  err = ncmpi_def_var(fileID, datasetName, ncDataType,
		      dims, dimIDs, &tmpVar);
  assert(err == NC_NOERR);
}
