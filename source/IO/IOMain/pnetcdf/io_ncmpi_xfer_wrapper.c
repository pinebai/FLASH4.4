#include "io_ncmpi_xfer_wrapper.h"

int io_ncmpi_xfer_wrapper(const int myPE, const int fileID, const int xferType,
			  const char datasetName[], const int memType,
			  const int memSize[], const int memStart[],
			  const int memCount[], const int diskStart[],
			  const int diskCount[], const int dims, void * pData)
{
  MPI_Offset mpiDiskStart[IO_MAX_DIMS], mpiDiskCount[IO_MAX_DIMS],
    mpiMemCountScalar;
  MPI_Datatype mpiMemTypeScalar, mpiMemType;
  int varID, err, i, ierr, usePrimitiveTypeXfer;
  int order = MPI_ORDER_C;

  err = ncmpi_inq_varid(fileID, datasetName, &varID);

  if (err != NC_NOERR) {

    if (xferType == IO_READ_XFER) {
      ierr = -1;
      if (myPE == MASTER_PE) {
	printf(" [%s]: Skipping missing variable '%s'.\n", __FILE__, datasetName);
      }
    } else {
      printf(" [%s]: Unable to locate variable %s.\n", __FILE__, datasetName);
      ierr = Driver_abortFlashC("Error! pnetcdf variable not defined");
    }

  } else {

    /* Describe the data in file */
    for (i=0; i<dims; ++i) {
      if (diskStart[i] > 0 && diskCount[i] == 0) {
        /* Setting mpiDiskStart[i] to 0 stops an overly cautious index out of
	   bounds assertion error in pnetcdf library >= version 1.2.0pre1. */
	mpiDiskStart[i] = 0;
      } else {
	mpiDiskStart[i] = (MPI_Offset) diskStart[i];
      }
      mpiDiskCount[i] = (MPI_Offset) diskCount[i];
    }


    /* Describe the data in memory - attempt to use primitive types first */
    usePrimitiveTypeXfer = 1;

    /* Look at all dimension except index 0 because that is the slowest
       varying dimension.  As long as memsize and memcount are the same for
       all except the slowest varying dimension we can use primitive types. */
    for (i=dims-1; i>0; --i) {
      if (memSize[i] != memCount[i]) {
	usePrimitiveTypeXfer = 0;
      }
    }

    /* memStart needs to be 0 because we will obtain the memory count using
       the expression: mpiMemCountScalar *= (MPI_Offset) memCount[i]; */
    for (i=0; i<dims; ++i) {
      if (diskCount[i] != memCount[i] || memStart[i] != 0) {
	usePrimitiveTypeXfer = 0;
      }
    }


    mpiMemTypeScalar = io_mpi_type_primitive(memType);
    mpiMemCountScalar = 1;

    if (0 == usePrimitiveTypeXfer) {
      /* Define a subarray to describe non-contiguous data in memory */
      if (myPE == MASTER_PE) {
	printf(" [%s]: Non-contiguous data in memory for dataset '%s'.\n",
	       __FILE__, datasetName);
      }
      err = MPI_Type_create_subarray(dims, (int*) memSize, (int*) memCount,
				     (int*) memStart, order, mpiMemTypeScalar,
				     &mpiMemType);
      assert(err == MPI_SUCCESS);
      err = MPI_Type_commit(&mpiMemType);
      assert(err == MPI_SUCCESS);
    } else {
      mpiMemType = mpiMemTypeScalar;
      for (i=0; i<dims; ++i) {
	mpiMemCountScalar *= (MPI_Offset) memCount[i];
      }
    }


    err = io_ncmpi_xfer(myPE, fileID, xferType, varID,
			mpiMemType, mpiMemCountScalar,
			mpiDiskStart, mpiDiskCount, dims, pData);
    assert (err == 0);


    if (0 == usePrimitiveTypeXfer) {
      err = MPI_Type_free(&mpiMemType);
      assert(err == MPI_SUCCESS);
    }

    ierr = 0;
  }

  return ierr;
}
