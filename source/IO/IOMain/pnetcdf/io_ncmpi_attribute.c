#include "io_ncmpi_attribute.h"

/* WARNING: The caller must null-terminate all strings.
   Must be in define mode before using either of these functions!!! */

void io_ncmpi_attribute_create(const int myPE,
			       const int fileID,
			       const int diskType,
			       const int dims,
			       const int datasetSize[],
			       const char datasetName[],
			       const char attDatasetName[])
{
  MPI_Offset datasetSizeScalar;
  int varID, err, i;
  nc_type ncDiskType;
  double *dblBuf = NULL;
  float *fltBuf = NULL;
  int *intBuf = NULL;
  char *chrBuf = NULL;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  assert(dims > 0);
  for (i=0; i<dims; ++i) {
    assert(datasetSize[i] > 0);
  }

  if (debugIO && myPE == MASTER_PE) {
    printf(" [%s]: Dataset: %s, creating attribute: %s.\n",
	   __FILE__, datasetName, attDatasetName);
  }

  if (strcmp(datasetName,"NC_GLOBAL") == 0) {
    varID = NC_GLOBAL;
  } else {
    err = ncmpi_inq_varid(fileID, datasetName, &varID);
    assert(err == NC_NOERR);
  }

  ncDiskType = io_ncmpi_type_primitive(diskType);
  datasetSizeScalar = 1;
  for (i=0; i<dims; ++i) {
    datasetSizeScalar *= datasetSize[i];
  }

  /* Allocate using calloc to stop false positives from valgrind */
  switch (ncDiskType) {
  case (NC_INT):
    intBuf = calloc (datasetSizeScalar, sizeof(*intBuf));
    err = ncmpi_put_att_int(fileID, varID, attDatasetName, ncDiskType,
			    datasetSizeScalar, intBuf);
    assert(err == NC_NOERR);
    free (intBuf);
    break;
  case (NC_DOUBLE):
    dblBuf = calloc (datasetSizeScalar, sizeof(*dblBuf));
    err = ncmpi_put_att_double(fileID, varID, attDatasetName, ncDiskType,
			       datasetSizeScalar, dblBuf);
    assert(err == NC_NOERR);
    free (dblBuf);
    break;
  case (NC_FLOAT):
    fltBuf = calloc (datasetSizeScalar, sizeof(*fltBuf));
    err = ncmpi_put_att_float(fileID, varID, attDatasetName, ncDiskType,
			      datasetSizeScalar, fltBuf);
    assert(err == NC_NOERR);
    free (fltBuf);
    break;
  case (NC_CHAR):
    chrBuf = calloc (datasetSizeScalar, sizeof(*chrBuf));
    err = ncmpi_put_att_text(fileID, varID, attDatasetName,
			     datasetSizeScalar, chrBuf);
    assert(err == NC_NOERR);
    free (chrBuf);
    break;
  default:
    Driver_abortFlashC("[io_ncmpi_attribute_create] Invalid file type");
  }
}


void io_ncmpi_attribute_write(const int myPE,
			      const int fileID,
			      const int memType,
			      const char datasetName[],
			      const char attDatasetName[],
			      const void * const pData)
{
  MPI_Offset memCount;
  int varID, err;
  nc_type ncDiskType, ncMemDatatype;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  if (debugIO && myPE == MASTER_PE) {
    printf(" [%s]: Dataset: %s, writing attribute: %s.\n",
	   __FILE__, datasetName, attDatasetName);
  }

  if (strcmp(datasetName,"NC_GLOBAL") == 0) {
    varID = NC_GLOBAL;
  } else {
    err = ncmpi_inq_varid(fileID, datasetName, &varID);
    assert(err == NC_NOERR);
  }

  err = ncmpi_inq_atttype(fileID, varID, attDatasetName, &ncDiskType);
  assert(err == NC_NOERR);

  err = ncmpi_inq_attlen(fileID, varID, attDatasetName, &memCount);
  assert(err == NC_NOERR);


  ncMemDatatype = io_ncmpi_type_primitive(memType);
  switch (ncMemDatatype) {
  case (NC_INT):
    err = ncmpi_put_att_int(fileID, varID, attDatasetName, ncDiskType,
			    memCount, pData);
    assert(err == NC_NOERR);
    break;
  case (NC_DOUBLE):
    err = ncmpi_put_att_double(fileID, varID, attDatasetName, ncDiskType,
			       memCount, pData);
    assert(err == NC_NOERR);
    break;
  case (NC_FLOAT):
    err = ncmpi_put_att_float(fileID, varID, attDatasetName, ncDiskType,
			      memCount, pData);
    assert(err == NC_NOERR);
    break;
  case (NC_CHAR):
    err = ncmpi_put_att_text(fileID, varID, attDatasetName,
			     memCount, pData);
    assert(err == NC_NOERR);
    break;
  default:
    Driver_abortFlashC("[io_ncmpi_attribute_write] Invalid memory type");
  }
}
