#include "Test_parallel_write.h"

#define NXBG (NXB+2*K1D*NGUARD)
#define NYBG (NYB+2*K2D*NGUARD)
#define NZBG (NZB+2*K3D*NGUARD)
#define VAR_LEN 4
#define VAR_LEN_PLUS_NULL (VAR_LEN + 1)

/* Use a 1D array with the following macro to simulate a 5D array */
#define SUBSCRIPT5D(b,z,y,x,u) ((b)*NZBG*NYBG*NXBG*NUNK_VARS) + \
  ((z)*NYBG*NXBG*NUNK_VARS) +					\
  ((y)*NXBG*NUNK_VARS) +					\
  ((x)*NUNK_VARS) +						\
  (u)

/* Test requires 4 processes with each process writing a single FLASH block */
#define NPROCS 4
#define NBLOCKS 1

int main(int argc, char *argv[])
{

  /* Metadata definitions:
     -------------------------------------------------------------------------*/
  const double minVar[NUNK_VARS] = {2.1, 5.3};
  const double maxVar[NUNK_VARS] = {3.2, 6.4};


  /* Data definitions:
     -------------------------------------------------------------------------*/
  double data[NBLOCKS * NZBG * NYBG * NXBG * NUNK_VARS];
  const int blockOuterSize[] = {NXBG, NYBG, NZBG};
  const int blockInnerSize[] = {NXB, NYB, NZB};
  const int blockInnerOffset[] = {NGUARD*K1D, NGUARD*K2D, NGUARD*K3D};
  const int numDataDims = IO_MESH_DIMS;
  int globalOffsetArray[IO_MESH_DIMS];
  int localSubSize[IO_MESH_DIMS];
  int attDims = 1;

  /* Other definitions:
     -------------------------------------------------------------------------*/
  const char filename[] = "Test_parallel_write.pnetcdf";
  const char *dsetNames[2] = {"UNK___Two_variables_in_double_precision",
			      "UNK___One_variable_in_single_precision"};
  const char *attNames[2] = {"minimum", "maximum"};
  const int attChkSize[1] = {2};
  const int attPltSize[1] = {1};
  const int dataType[2] = {IO_FLASH_DOUBLE, IO_FLASH_FLOAT};
  const int memType[2] = {IO_FLASH_DOUBLE, IO_FLASH_DOUBLE};
  const int gridDataStructs[5] = {CENTER, FACEX, FACEY, FACEZ, SCRATCH};
  int gridDataStruct;
  int numVar, numPlotVar;
  const int plotVarArr[1] = {0};
  int fileID, mpiRank, mpiSize, i, b, z, y, x, u, v;
  int dim_tot_blocks, dim_nxb, dim_nyb, dim_nzb, dim_nvar, err, cumulative;
  int dimids[IO_MESH_DIMS], dim_onevar, fileType;


  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  if (mpiSize != NPROCS) {
    Driver_abortFlashC("This is a 4-processor unit test");
  }
  printf("Hello from processor %d\n", mpiRank);


  /* Data initialisation:
     -------------------------------------------------------------------------*/
  /* Fill with garbage */
  cumulative = 0;
  for (b=0; b<NBLOCKS; ++b) {
    for (z=0; z<(NZB+2*NGUARD*K3D); ++z) {
      for (y=0; y<(NYB+2*NGUARD*K2D); ++y) {
	for (x=0; x<(NXB+2*NGUARD*K1D); ++x) {
	  for (u=0; u<NUNK_VARS; ++u) {
	    data [ SUBSCRIPT5D(b,z,y,x,u) ] = -(double) cumulative;
	    cumulative = cumulative + 1;
	  }
	}
      }
    }
  }


  /* Fill internal cells with mpiRank value */
  for (b=0; b<NBLOCKS; ++b) {
    for (z=NGUARD*K3D; z<(NZB+NGUARD*K3D); ++z) {
      for (y=NGUARD*K2D; y<(NYB+NGUARD*K2D); ++y) {
	for (x=NGUARD*K1D; x<(NXB+NGUARD*K1D); ++x) {
	  for (u=0; u<NUNK_VARS; ++u) {
	    data [ SUBSCRIPT5D(b,z,y,x,u) ] = (double) mpiRank;
	  }
	}
      }
    }
  }



  /* Create a Pnetcdf file */
  err = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER,
		     MPI_INFO_NULL, &fileID);
  assert(err == NC_NOERR);


  /* Pnetcdf define mode */
  /* ------------------------------------------------------------------------ */
  err = ncmpi_def_dim(fileID, "dim_tot_blocks", (MPI_Offset)(NBLOCKS*NPROCS), &dim_tot_blocks);
  assert(err == NC_NOERR);
  err = ncmpi_def_dim(fileID, "dim_nzb", (MPI_Offset)(NZB), &dim_nzb);
  assert(err == NC_NOERR);
  err = ncmpi_def_dim(fileID, "dim_nyb", (MPI_Offset)(NYB), &dim_nyb);
  assert(err == NC_NOERR);
  err = ncmpi_def_dim(fileID, "dim_nxb", (MPI_Offset)(NXB), &dim_nxb);
  assert(err == NC_NOERR);
  err = ncmpi_def_dim(fileID, "dim_nvar", (MPI_Offset)(NUNK_VARS), &dim_nvar);
  assert(err == NC_NOERR);
  err = ncmpi_def_dim(fileID, "dim_onevar", (MPI_Offset)(1), &dim_onevar);
  assert(err == NC_NOERR);


  /* Global disk size */
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_nzb;
  dimids[2] = dim_nyb;
  dimids[3] = dim_nxb;
  dimids[4] = dim_nvar;
  io_ncmpi_create_dataset(mpiRank,
			  fileID,
			  dataType[0],
			  numDataDims,
			  dimids,
			  dsetNames[0]);

  io_ncmpi_attribute_create(mpiRank,
			    fileID,
			    dataType[0],
			    attDims,
			    attChkSize,
			    dsetNames[0],
			    attNames[0]);

  io_ncmpi_attribute_create(mpiRank,
			    fileID,
			    dataType[0],
			    attDims,
			    attChkSize,
			    dsetNames[0],
			    attNames[1]);


  dimids[4] = dim_onevar;
  io_ncmpi_create_dataset(mpiRank,
			  fileID,
			  dataType[1],
			  numDataDims,
			  dimids,
			  dsetNames[1]);

  io_ncmpi_attribute_create(mpiRank,
			    fileID,
			    dataType[1],
			    attDims,
			    attPltSize,
			    dsetNames[1],
			    attNames[0]);

  io_ncmpi_attribute_create(mpiRank,
			    fileID,
			    dataType[1],
			    attDims,
			    attPltSize,
			    dsetNames[1],
			    attNames[1]);


  err = ncmpi_enddef(fileID);
  assert(err == NC_NOERR);
  /* ------------------------------------------------------------------------ */



  /* After all this setup we use the following 3 functions to write
     the data and amended attributes to disk */
  /* ------------------------------------------------------------------------ */
  numVar = NUNK_VARS;
  for (i=0; i<5; ++i) {
    gridDataStruct = gridDataStructs[i];
    if (gridDataStruct == CENTER) {
      numVar = NUNK_VARS;
      numPlotVar = 1;
    } else {
      numVar = 0;
      numPlotVar = 0;
    }
    /* Initialize all types even though we will just be using CENTER - requried
       because the free routine releases memory for all data structure */
    FTOC(io_init_grid_mpi_types)(&mpiRank,
				 &gridDataStruct,
				 &blockOuterSize[0],
				 &blockInnerSize[0],
				 &blockInnerOffset[0],
				 &numVar,
				 &plotVarArr[0],
				 &numPlotVar);
  }

  numVar = NUNK_VARS;
  numPlotVar = 1;

  /* Total amount of local data to write: */
  localSubSize[0] = NBLOCKS;
  localSubSize[1] = blockInnerSize[2];
  localSubSize[2] = blockInnerSize[1];
  localSubSize[3] = blockInnerSize[0];
  localSubSize[4] = numVar;

  /* MyPE's global offset in the data file: */
  globalOffsetArray[0] = mpiRank;
  globalOffsetArray[1] = 0;
  globalOffsetArray[2] = 0;
  globalOffsetArray[3] = 0;
  globalOffsetArray[4] = 0;


  gridDataStruct = CENTER;
  for (v=0; v<2; ++v) {
    if (v == 0) {
      fileType = CHECKPOINTFILE;
      numVar = NUNK_VARS;
    }
    if (v == 1) {
      fileType = PLOTFILE;
      numVar = numPlotVar;
    }

    localSubSize[4] = numVar;


    /* Call the actual function that writes the data */
    io_ncmpi_xfer_mesh_dataset(mpiRank,
			       fileID,
			       IO_WRITE_XFER,
			       fileType,
			       10,
			       gridDataStruct,
			       numDataDims,
			       0,
			       globalOffsetArray,
			       localSubSize,
			       dsetNames[v],
			       &data[0]);

    io_ncmpi_attribute_write(mpiRank,
			     fileID,
			     memType[v],
			     dsetNames[v],
			     attNames[0],
			     minVar);

    io_ncmpi_attribute_write(mpiRank,
			     fileID,
			     memType[v],
			     dsetNames[v],
			     attNames[1],
			     maxVar);
  }


  FTOC(io_free_grid_mpi_types)();
  /* ------------------------------------------------------------------------ */

  err = ncmpi_close(fileID);
  assert(err == NC_NOERR);

  MPI_Finalize();
  return 0;
}
