#include "io_set_grid_mpi_types.h"

#ifdef INDEXREORDER
# error [FLASH message] Derived datatype I/O incompatible with -index-reorder!
#endif

MPI_Datatype G_mesh_all_var_type[NUM_MESH_TYPES],
  G_mesh_subset_var_type[NUM_MESH_TYPES],
  G_mesh_single_var_type[NUM_MESH_TYPES];

void FTOC(io_init_grid_mpi_types)(const int * const pMyPE,
				  const int * const pGridDataStruct,
				  const int blockOuterSize[],
				  const int blockInnerSize[],
				  const int blockInnerOffset[],
				  const int * const pNumGridVar,
				  const int plotVarArr[],
				  const int * const pNumPlotVar)
{
  const int myPE = *pMyPE;
  const int gridDataStruct = *pGridDataStruct;
  const int numGridVar = *pNumGridVar;
  const int numPlotVar = *pNumPlotVar;
  int err, i, count, index, sizeDouble;

  MPI_Aint printNewExtent;
  int ndims = IO_MESH_DIMS-1;
  int array_of_sizes[IO_MESH_DIMS-1], array_of_subsizes[IO_MESH_DIMS-1],
    array_of_starts[IO_MESH_DIMS-1], array_of_single_var_subsizes[IO_MESH_DIMS-1];
  int order = MPI_ORDER_C;

  const char * struct_names[IO_MESH_DIMS+1] =
    {"NULL","scratch","unk","facex","facey","facez"};
  int ndims_plot = ndims-1; /* Use 3 so we exclude NUNK_VARS dimension */
  int array_of_blocklengths[MAX_MESH_VAR], array_of_displacements[MAX_MESH_VAR];
  MPI_Datatype varPickup, varStencil, *pAllVarType = NULL,
    *pSubsetVarType = NULL, *pSingleVarType = NULL;
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  assert(numPlotVar >= 0 && numPlotVar <= MAX_MESH_VAR);
  assert(numGridVar >= 0 && numGridVar <= MAX_MESH_VAR);
  assert(ndims == 4);
  assert(ndims_plot == 3);

  /* These values are required so that the strings in struct_names
     correspond to the integer grid structure values */
  assert(SCRATCH == 1);
  assert(CENTER == 2);
  assert(FACEX == 3);
  assert(FACEY == 4);
  assert(FACEZ == 5);
  assert(gridDataStruct >= 1 && gridDataStruct <= 5);

  pAllVarType = &G_mesh_all_var_type[gridDataStruct-1];
  pSubsetVarType = &G_mesh_subset_var_type[gridDataStruct-1];
  pSingleVarType = &G_mesh_single_var_type[gridDataStruct-1];


  /* The final datatypes will describe the internal region of a block.
     Copy block information into arrays suitable for MPI_Type_create_subarray.
     The block interior for checkpoint files is described by 3 4-dimensional
     arrays.  However, we can re-use these arrays for plotfiles by just
     considering the first 3 dimensions. */

  /* The comments after each line indicate the pre-processor values for UNK */
  array_of_sizes[0] = blockOuterSize[2]; /* NZB+2*GZ */
  array_of_sizes[1] = blockOuterSize[1]; /* NYB+2*GY */
  array_of_sizes[2] = blockOuterSize[0]; /* NXB+2*GX */
  array_of_sizes[3] = numGridVar; /* NUNK_VARS */

  array_of_subsizes[0] = blockInnerSize[2]; /* NZB */
  array_of_subsizes[1] = blockInnerSize[1]; /* NYB */
  array_of_subsizes[2] = blockInnerSize[0]; /* NXB */
  array_of_subsizes[3] = numGridVar; /* NUNK_VARS */

  array_of_single_var_subsizes[0] = array_of_subsizes[0]; /* NZB */
  array_of_single_var_subsizes[1] = array_of_subsizes[1]; /* NYB */
  array_of_single_var_subsizes[2] = array_of_subsizes[2]; /* NXB */
  array_of_single_var_subsizes[3] = 1; /* 1 */

  array_of_starts[0] = blockInnerOffset[2]; /* NGUARD*K3D-1 */
  array_of_starts[1] = blockInnerOffset[1]; /* NGUARD*K2D-1 */
  array_of_starts[2] = blockInnerOffset[0]; /* NGUARD-1 */
  array_of_starts[3] = 0; /* 0 */


  /*
     1. Create a type describing all internal variables.
     -----------------------------------------------------------------------
  */
  if (numGridVar > 0) {
    /* Define a subarray to exclude guardcells */
    err = MPI_Type_create_subarray(ndims, array_of_sizes, array_of_subsizes,
				   array_of_starts, order, MPI_DOUBLE_PRECISION,
				   pAllVarType);
    assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(pAllVarType);
    assert(err == MPI_SUCCESS);
  } else {
    *pAllVarType = MPI_DATATYPE_NULL;
  }
  /* ----------------------------------------------------------------------- */


  /*
     2. Create a type describing a single internal variable.
     -----------------------------------------------------------------------
  */
  if (numGridVar > 0) {
    err = MPI_Type_create_subarray(ndims, array_of_sizes,
				   array_of_single_var_subsizes,
				   array_of_starts, order, MPI_DOUBLE_PRECISION,
				   pSingleVarType);
    assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(pSingleVarType);
    assert(err == MPI_SUCCESS);
  } else {
    *pSingleVarType = MPI_DATATYPE_NULL;
  }
  /* ----------------------------------------------------------------------- */


  /*
     3. Create a type describing a subset of internal variables.
     -----------------------------------------------------------------------
  */
  if (debugIO) {
    if (myPE == MASTER_PE) {
      printf(" [%s]: Mesh structure %s - describing %d plot variables"
	     " (selected from %d mesh variables).\n",
	     __FILE__, struct_names[gridDataStruct], numPlotVar, numGridVar);
    }
  }


  /* Important! We assume plotVarArr is sorted and holds zero-based data */
  if (numPlotVar > 0) {
    index = 0;
    array_of_displacements[0] = plotVarArr[0];
    array_of_blocklengths[0] = 1;

    /* If there is more than 1 contiguous variable we try to
       combine the variables into larger blocks. */
    if (numPlotVar > 1) {

      for (i=1; i<numPlotVar; ++i) {

	/* Test for two variables being next to each other */
	if ( (array_of_displacements[index] + array_of_blocklengths[index])
	     == plotVarArr[i]) {

	  array_of_blocklengths[index] = array_of_blocklengths[index] + 1;

	} else {

	  /* These variables are not next to each other so start a new block */
	  index = index + 1;
	  array_of_displacements[index] = plotVarArr[i];
	  array_of_blocklengths[index] = 1;
	}
      }
    }
    count = index + 1;

    if (debugIO) {
      if (myPE == MASTER_PE) {
	for (i=0; i<count; ++i) {
	  printf(" [%s]: array_of_displacements[%d]=%d, "
		 "array_of_blocklengths[%d]=%d\n", __FILE__,
		 i, array_of_displacements[i], i, array_of_blocklengths[i]);
	}
      }
    }


    /* First create a type to pick up the data we need */
    err = MPI_Type_indexed(count, array_of_blocklengths,
			   array_of_displacements, MPI_DOUBLE_PRECISION, &varPickup);
    assert(err == MPI_SUCCESS);

    /* Now adjust the extents */
    err = MPI_Type_size(MPI_DOUBLE_PRECISION, &sizeDouble);
    assert(err == MPI_SUCCESS);


    /* This is a simple wrapper around MPI_Type_create_resized unless
       we are using an MPI-1 MPI implementation. */
    err = io_resize_mpi_type(varPickup, (MPI_Aint)0,
			     (MPI_Aint)(numGridVar * sizeDouble), &varStencil);
    assert(err == MPI_SUCCESS);

    if (debugIO) {
      if (myPE == MASTER_PE) {
	err = MPI_Type_extent(varStencil,&printNewExtent);
	assert(err == MPI_SUCCESS);
	printf(" [%s]: A single mesh variable is %d bytes."
	       " Extent of MPI type to select plot variables is %d bytes.\n",
	       __FILE__, sizeDouble, (int)printNewExtent);
      }
    }


    /* Now define a subarray to exclude guardcells */
    err = MPI_Type_create_subarray(ndims_plot, array_of_sizes, array_of_subsizes,
				   array_of_starts, order, varStencil,
				   pSubsetVarType);
    assert(err == MPI_SUCCESS);
    err = MPI_Type_commit(pSubsetVarType);
    assert(err == MPI_SUCCESS);

  } else {

    /* numPlotVar is zero.  Set global datatype to null datatype so the
       application crashes early if we mistakenly use the datatype. */
    *pSubsetVarType = MPI_DATATYPE_NULL;

  }
  /* ----------------------------------------------------------------------- */


}


void FTOC(io_free_grid_mpi_types)(void)
{
  int err, i;

  /* Free all derived data types. */
  for (i=0; i<NUM_MESH_TYPES; ++i) {
    if (G_mesh_all_var_type[i] != MPI_DATATYPE_NULL) {
      err = MPI_Type_free(&G_mesh_all_var_type[i]);
      assert(err == MPI_SUCCESS);
    }

    if (G_mesh_subset_var_type[i] != MPI_DATATYPE_NULL) {
      err = MPI_Type_free(&G_mesh_subset_var_type[i]);
      assert(err == MPI_SUCCESS);
    }

    if (G_mesh_single_var_type[i] != MPI_DATATYPE_NULL) {
      err = MPI_Type_free(&G_mesh_single_var_type[i]);
      assert(err == MPI_SUCCESS);
    }
  }
}


int io_resize_mpi_type(MPI_Datatype oldtype,
		       MPI_Aint lb,
		       MPI_Aint extent,
		       MPI_Datatype *pNewtype)
{
  int err;

#if MPI_VERSION >= 2

  /* This version uses the right tool for the job.  It uses a newer
     MPI function that is only available in MPI-2 implementations */
  err = MPI_Type_create_resized(oldtype, lb, extent, pNewtype);
  assert(err == MPI_SUCCESS);

#else

  /* This version uses the deprecated MPI_UB marker method to adjust the
     extent of an MPI datatype.  It is provided becase MPI_Type_create_resized
     is not available in MPI-1 implementations such as mpich-1.2.7p1.  The
     MPI_UB method is not a good solution compared to MPI_Type_create_resized
     and so it is only used as a last resort */

  MPI_Datatype array_of_types[2];
  MPI_Aint array_of_displacements[2];
  int array_of_blocklens[2];
  int count;

  /* I do not support modification of MPI_LB of the new type */
  assert(lb == (MPI_Aint)0);

  count = 2;
  array_of_blocklens[0] = 1;  /* 1 element of oldtype */
  array_of_blocklens[1] = 1;  /* 1 element of MPI_UB */

  array_of_displacements[0] = 0;
  array_of_displacements[1] = extent;

  array_of_types[0] = oldtype;
  array_of_types[1] = MPI_UB;

  err = MPI_Type_struct(count, array_of_blocklens, array_of_displacements,
			array_of_types, pNewtype);
  assert(err == MPI_SUCCESS);

#endif

  return err;
}
