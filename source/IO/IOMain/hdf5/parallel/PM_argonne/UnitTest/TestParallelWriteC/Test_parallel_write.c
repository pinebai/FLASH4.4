#include "Test_parallel_write.h"

/* We will have this many variables. */
#define NUNK_VARS 4

/* We will have this many cells (= numInternalCells+(2*numGuardCells)) */
#define NUM_CELLS 8

#define NUM_FILES 2

int HDF5_MODE = 0;

int main(int argc, char *argv[])
{

  /* Metadata definitions:
     -------------------------------------------------------------------------*/
  const double minVar[NUNK_VARS] =
    {2.1, 5.3, 8.5, 9.7};
  const double maxVar[NUNK_VARS] =
    {3.2, 6.4, 9.6, 10.8};
  int numVar;

  /* Data definitions:
     -------------------------------------------------------------------------*/
  const int numInternalCells = 4;
  const int numGuardCells = 2;
  const int numDataDims = 1;
  double data[NUM_CELLS] =
    {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
  const char attMinName[] = "minimum";
  const char attMaxName[] = "maximum";

  /* Other definitions:
     -------------------------------------------------------------------------*/
  const int useCollectiveHDF5 = 1;
  const char * filenames[NUM_FILES] = {"CheckpointFile.hdf5", "PlotFile.hdf5"};
  const char * datasetName = "unknown";
  const char * filename;
  const int diskType[NUM_FILES] = {IO_FLASH_DOUBLE, IO_FLASH_FLOAT};
  const int memType = IO_FLASH_DOUBLE;
  const int xferType = IO_WRITE_XFER;
  const int typeMatchedXfer = 1;

  hid_t fileAccessListID;
  herr_t err;
  int fileID; /* Should be hid_t, but will be int from calling Fortran code */
  int mpiRank, mpiSize, i, f, fileType;
  int memSize[1], memStart[1], memCount[1], diskSize[1], diskStart[1],
    diskCount[1];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if (mpiSize != 4) {
    Driver_abortFlashC("This is a 4-processor unit test");
  }
  printf("Hello from processor %d\n", mpiRank);


  HDF5_MODE = useCollectiveHDF5;


  /* Data initialisation:
     Each processor has a 1D-array consisting of 8 elements (4 internal
     elements and 2 guards on left and 2 guard cells on right).
     The unit-test runs on 4 processors and outputs
     a 16-element global 1D-array: 0000 1111 2222 3333
     -------------------------------------------------------------------------*/
  for (i=numGuardCells; i<numInternalCells+numGuardCells; ++i) {
    data[i] = mpiRank;
  }


  for (f=0; f<NUM_FILES; ++f) {

    filename = filenames[f];
    if (f == 0) {
      fileType = CHECKPOINTFILE;
      numVar = 4;
      if (mpiRank == MASTER_PE) {
	printf("We expect checkpoint file to contain: 0000 1111 2222 3333.\n");
      }
    }
    if (f == 1) {
      fileType = PLOTFILE;
      numVar = 2;
      if (mpiRank == MASTER_PE) {
	printf("We expect plot file to contain: 00 11 22 33.\n");
      }
    }


    /* PARALLEL IO DETAILS COPIED FROM:
       http://www.hdfgroup.org/HDF5/Tutor/examples/parallel/
       ------------------------------------------------------------------------ */
    /* Set up file access property list with parallel I/O access */
    fileAccessListID = H5Pcreate(H5P_FILE_ACCESS);
    assert(fileAccessListID != -1);
    err = H5Pset_fapl_mpio(fileAccessListID, MPI_COMM_WORLD, MPI_INFO_NULL);
    assert(err >= 0);

    /* Create a new file collectively and release property list identifier. */
    fileID = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fileAccessListID);
    assert(fileID >= 0);

    err = H5Pclose(fileAccessListID);
    assert(err >= 0);
  /* ------------------------------------------------------------------------ */


    memSize[0] = NUM_CELLS;
    memStart[0] = numGuardCells;
    memCount[0] = numVar;
    diskSize[0] = mpiSize * numVar; /* {16} or {8} */
    diskStart[0] = mpiRank * numVar; /* {0, 4, 8, 12} or {0, 2, 4, 6} */
    diskCount[0] = numVar;


    /* Create the single dataset and the 2 attributes */
    io_h5create_dataset(mpiRank,
			fileID,
			diskType[f],
			numDataDims,
			diskSize,
			datasetName);

    io_h5_attribute_create(mpiRank,
			   fileID,
			   diskType[f],
			   numDataDims,
			   diskCount,
			   datasetName,
			   attMinName);

    io_h5_attribute_create(mpiRank,
			   fileID,
			   diskType[f],
			   numDataDims,
			   diskCount,
			   datasetName,
			   attMaxName);


    /* Write the data and attributes */
    io_h5_xfer_wrapper(mpiRank,
		       fileID,
		       xferType,
		       typeMatchedXfer,
		       datasetName,
		       memType,
		       memSize,
		       memStart,
		       memCount,
		       diskStart,
		       diskCount,
		       numDataDims,
		       data);

    io_h5_attribute_write(mpiRank,
			  fileID,
			  memType,
			  datasetName,
			  attMinName,
			  minVar);

    io_h5_attribute_write(mpiRank,
			  fileID,
			  memType,
			  datasetName,
			  attMaxName,
			  maxVar);


    err = H5Fclose(fileID);
    assert(err >= 0);

  } /* End of for-each file */

  MPI_Finalize();
  return 0;
}
