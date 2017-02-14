#include "io_h5_xfer_mesh_dataset.h"

void io_h5_xfer_mesh_dataset(const int myPE,
			     const int fileID,
			     const int xferType,
			     const int fileFmt,
			     const int fileType,
			     const int gridDataStruct,
			     const int numFileDims,
			     const int numGridVar,
			     const int outputVarOffset[],
			     const int globalOffset[],
			     const int localSize[],
			     const int localSubSize[],
			     const int localOffset[],
			     const char datasetName[],
			     double * pData)
{
  const hid_t hFileID = (hid_t) fileID;
  herr_t (*old_func)(void*) = NULL;
  void *old_client_data = NULL;
  hsize_t memStart[IO_MESH_DIMS], diskStart[IO_MESH_DIMS],
    memDimens[IO_MESH_DIMS],
    memCount[IO_MESH_DIMS], diskCount[IO_MESH_DIMS];
  hsize_t memDummyDimens[IO_MESH_DIMS] = {1,1,1,1,1};
  hid_t hXferList, dataSpace, memSpace, dataSet;
  herr_t err;
  size_t datasetNameLen;
  int i, j, singleHyperslab, ierr, isNullSelection;
  char *fileTypeText = NULL, *xferTypeText = NULL;
  char legacyDatasetName[LEGACY_LABEL_LEN+1];
  const char *spacePad = "    "; /* 4 characters */
#ifdef DEBUG_IO
  const int debugIO = 1;
#else
  const int debugIO = 0;
#endif

  /* Check input variables are sensible */
  assert(numFileDims <= IO_MESH_DIMS);
  assert(fileFmt <= 10); /* Invalid file fmt */
  assert(numGridVar >= 1 && numGridVar <= MAX_MESH_VAR); /* Invalid num vars */
  assert(myPE >= 0); /* Invalid process identifier */
  assert(fileType == CHECKPOINTFILE || fileType == PLOTFILE);
  assert(gridDataStruct == CENTER || gridDataStruct == FACEX ||
	 gridDataStruct == FACEY || gridDataStruct == FACEZ ||
	 gridDataStruct == SCRATCH);

  if (debugIO) {
    if (fileType == CHECKPOINTFILE) fileTypeText = "checkpoint file";
    if (fileType == PLOTFILE) fileTypeText = "plot file";
    if (xferType == IO_WRITE_XFER) xferTypeText = "write";
    if (xferType == IO_READ_XFER) xferTypeText = "read";
    if (myPE == MASTER_PE) {
      printf(" [%s]: About to %s dataset %s in %s.\n",
	     __FILE__, xferTypeText, datasetName, fileTypeText);
    }
  }

  /* open the dataset
     suppress the error message if the variable does not exist in file */
  err = H5Eget_auto(&old_func, &old_client_data);
  assert(err >= 0);
  err = H5Eset_auto(NULL, NULL);
  assert(err >= 0);

  dataSet = H5Dopen(hFileID, datasetName);
  /* If the open fails and this is a restart then we attempt to open
     the dataset again with a backwards-compatible space padded name.
     We do this because the original I/O in FLASH creates datasets
     which have names of 4 characters in length, even if the space
     trimmed name is less than 4 characters, e.g. 'al' is 'al  '  */
  if (dataSet < 0 && xferType == IO_READ_XFER) {
    datasetNameLen = strlen(datasetName);
    if (datasetNameLen < LEGACY_LABEL_LEN) {
      /* The if statement above ensures that the strncpy call below
	 will copy a null character into legacyDatasetName */
      strncpy(legacyDatasetName, datasetName, LEGACY_LABEL_LEN);
      strncat(legacyDatasetName, spacePad, LEGACY_LABEL_LEN - datasetNameLen);
      dataSet = H5Dopen(hFileID, legacyDatasetName);
    }
  }

  err = H5Eset_auto(old_func, old_client_data);
  assert(err >= 0);


  if (dataSet < 0) {

    if (xferType == IO_READ_XFER) {
      if (myPE == MASTER_PE) {
	printf(" [%s]: Couldn't find variable '%s' in file, so skipping it.\n",
	       __FILE__, datasetName);
      }
    } else {
      ierr = Driver_abortFlashC("Error! H5Dopen failed");
    }

  } else {

    dataSpace = H5Dget_space(dataSet);
    assert(dataSpace >= 0);

    isNullSelection = 0;
    for (i=0; i<IO_MESH_DIMS; ++i) {
      isNullSelection = isNullSelection || (localSubSize[i] == 0);
    }

    if (isNullSelection) {
      /* There are zero blocks or zero variables so we make a null selection */
      memSpace = H5Screate_simple(IO_MESH_DIMS, memDummyDimens, NULL);
      assert(memSpace >= 0);
      err = H5Sselect_none(memSpace);
      assert (err >= 0);
      err = H5Sselect_none(dataSpace);
      assert (err >= 0);
    } else {

      /* Construct arrays that descibe how data is stored locally and
	 where we want to store the data in the file. */
      for (i=0; i<IO_MESH_DIMS; ++i) {
	memDimens[i] = localSize[i];
	memStart[i] = localOffset[i]; /* Local offset in memory */
	memCount[i] = localSubSize[i];
      }

      for (i=0; i<numFileDims; ++i) {
	diskStart[i] = globalOffset[i]; /* Global offset on disk */
	diskCount[i] = localSubSize[i];
      }

      /* Define memSpace (how stored in memory) and hyperslab in memory.
	 The hyperslab selects the internal cells of a block (i.e. minus guard
	 cells).  In the case of plot files and file format 10 we use a
	 hyperslab to select the internal cells of a block for a subset
	 of FLASH variables. */
      memSpace = H5Screate_simple(IO_MESH_DIMS, memDimens, NULL);
      assert(memSpace >= 0);

      singleHyperslab =
	((fileFmt <= 9) ||
	 (fileFmt == 10 && fileType == CHECKPOINTFILE));

      if (singleHyperslab) {

	if (fileFmt <= 9) {
	  /* The address of the mesh pointer is already on the variable
	     we wish to output, and so memStart is 0 and not
	     localOffset[IO_MESH_DIMS-1] */
	  memStart[IO_MESH_DIMS-1] = 0;
	  memCount[IO_MESH_DIMS-1] = 1;  /* 1 mesh variable */
	}

	err = H5Sselect_hyperslab(memSpace, H5S_SELECT_SET,
				  memStart, NULL, memCount, NULL);
	assert(err >= 0);

      } else {

	/* Dump all internal data grid data for a subset of variables.  Strategy
	   is to create a hyperslab over one variable at a time, and if there
	   is more than 1 variable, successively add the hyperslabs. */
	memCount[IO_MESH_DIMS-1] = 1; /* Update "count" to contain 1 grid variable */
	err = H5Sselect_none(memSpace);
	assert(err >= 0);

	for (j=0; j<numGridVar; ++j) {
	  /* The address of the mesh pointer is on variable 0, and so we need
	     to use localOffset to reach variable N. */
	  memStart[IO_MESH_DIMS-1] = localOffset[IO_MESH_DIMS-1] + outputVarOffset[j];

	  /* Add this grid variable to the memory datatype */
	  err = H5Sselect_hyperslab(memSpace, H5S_SELECT_OR,
				    memStart, NULL, memCount, NULL);
	  assert(err >= 0);
	}
      }

      /* Create a hyperslab to describe myPE's portion in the file */
      err = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET,
				diskStart, NULL, diskCount, NULL);
      assert(err >= 0);
    }


    /* Set independent/collective transfer mode */
    hXferList = H5Pcreate(H5P_DATASET_XFER);
    assert (hXferList != -1);
#ifdef H5_HAVE_PARALLEL
    if (HDF5_MODE == COLLECTIVE) {
      err = H5Pset_dxpl_mpio(hXferList, H5FD_MPIO_COLLECTIVE);
      assert (err >= 0);
    } else {
      err = H5Pset_dxpl_mpio(hXferList, H5FD_MPIO_INDEPENDENT);
      assert (err >= 0);
    }
#endif



    if (xferType == IO_READ_XFER) {
      /* Read data from disk hyperslab to memory */
      err = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
		    hXferList, pData);
    } else if (xferType == IO_WRITE_XFER) {
      /* Write data from memory to disk hyperslab */
      err = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, memSpace, dataSpace,
		     hXferList, pData);
    }
    assert (err >= 0);

    if (debugIO) {
      err = io_h5_report_xfer_method(myPE, hXferList, datasetName);
      assert (err == 0);
    }


    err = H5Pclose(hXferList);
    assert (err >= 0);

    err = H5Sclose(memSpace);
    assert(err >= 0);

    err = H5Sclose(dataSpace);
    assert(err >= 0);

    err = H5Dclose(dataSet);
    assert(err >= 0);
  }
}
