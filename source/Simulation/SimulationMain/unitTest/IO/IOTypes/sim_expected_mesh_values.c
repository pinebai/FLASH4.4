#include "sim_expected_mesh_values.h"

/* After writing the output files, verify that the checkpoint file
   contains expected values. */

int FTOC(sim_expected_mesh_values)(const int * const pMyPE,
				   const char * const pLabels)
{
  const int myPE = *pMyPE;

  double att[NUNK_VARS], *buf = NULL;
  obj_size_t bufsize, bufptr, dims[5]; /* Up to 5D dataset */
  obj_handle_t fileID, dsetID;
  int i, a, b, offset, fail, nReads, nVar, fileFmt, ndims;
  char dsetStr[MAX_STRING_LENGTH];
  const char ExpDatasetName[] = "unknown";
  const char * AttDatasetName[NUM_ATTR] = {"minimum","maximum"};

#if defined (EXP_HDF5)
  const char FileName[] = "flash_hdf5_chk_0001";
#elif defined (EXP_PNETCDF)
  const char FileName[] = "flash_ncmpi_chk_0001";
#else
  const char FileName[] = "error";
  if (myPE == 0) {
    printf("Comparison test program designed for type based I/O only\n");
  }
  fail = SIM_FAIL;
  goto error;
#endif


  fail = SIM_SUCCESS;
  if (myPE == 0) {
    fileID = Obtain_File_Handle(FileName);
 
    fileFmt = Query_File_Fmt(myPE, fileID);
    nReads = (fileFmt == 9 ? NUNK_VARS : 1);
    nVar = (fileFmt == 9 ? 1 : NUNK_VARS);

    offset = 0;
    for (i=0; i<nReads; ++i) {

      /* Construct dataset name then open dataset */
      if (fileFmt == 9) {      
	memcpy(dsetStr, pLabels+offset, LABEL_LEN);
	dsetStr[LABEL_LEN] = '\0';
	offset += LABEL_LEN;
      } else {
	strcpy(dsetStr, ExpDatasetName);
      }
      dsetID = Obtain_Dataset_Handle(fileID, dsetStr);


      /* Find the extent of the dataset */
      Query_Dataset_extent(fileID, dsetID, fileFmt, dims, &ndims);
      bufsize = 1;
      for (b=0; b<ndims; ++b) {
	bufsize*= dims[b];
      }
      assert(bufsize > 0);


      /* Read dataset and check we get expected values */
      buf = malloc ((size_t)bufsize * sizeof(*buf));
      assert(buf != NULL);
      Read_Dataset(fileID, dsetID, dims, bufsize, buf);
      bufptr = 0;
      while (bufptr < bufsize) {
	if ((fail = check_mesh_variables(&buf[bufptr], i, nVar)) == SIM_FAIL) {
	  goto error;
	}
	bufptr += nVar;
      }
      free(buf);
      buf = NULL;


      /* Read "minimum" and then "maximum" attribute */
      for (a=0; a<NUM_ATTR; ++a) {
	Read_Attribute(fileID, dsetID, AttDatasetName[a], att);
	/* Check attribute contains expected value(s) */
	if ((fail = check_mesh_variables(&att[0], i, nVar)) == SIM_FAIL)
	  goto error;
      } /* End of NUM_ATTR attributes */

      Release_Dataset_Handle(dsetID);
    } /* End of nReads of dataset(s) */

    Release_File_Handle(fileID);
  } /* End of myPE == 0 */


error:
  /* Should really free resources, but this is not neccessary for such
     a simple piece of code.  I will not call I/O libraries anymore 
     and the program will terminate pretty much as soon as we return,
     so it will not cause a problem. */
  return fail;
}


obj_handle_t Obtain_File_Handle(const char FileName[])
{
  obj_handle_t fileID;
  obj_err_t err;

#if defined (EXP_HDF5)
  fileID = H5Fopen(FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(fileID >= 0);
#elif defined (EXP_PNETCDF)
  err = ncmpi_open(MPI_COMM_SELF, FileName, NC_NOWRITE, 
                   MPI_INFO_NULL, &fileID);
  assert(err == NC_NOERR);
#endif

  return fileID;
}


void Release_File_Handle(obj_handle_t fileID)
{
  obj_err_t err;

#if defined (EXP_HDF5)
  err = H5Fclose(fileID);
  assert(err >= 0);
#elif defined (EXP_PNETCDF)
  err = ncmpi_close(fileID);
  assert(err == NC_NOERR);
#endif
}


int Query_File_Fmt(const int myPE,
		   const obj_handle_t fileID)
{
  int fileFmt;

#if defined (EXP_HDF5)
  FTOC(io_h5_read_file_format)(&myPE, &fileID, &fileFmt);
#elif defined (EXP_PNETCDF)
  FTOC(io_ncmpi_read_file_format)(&myPE, &fileID, &fileFmt);
#endif
  assert(fileFmt == 9 || fileFmt == 10);

  return fileFmt;
}


obj_handle_t Obtain_Dataset_Handle(const obj_handle_t fileID,
				   const char dsetStr[])
{
  obj_handle_t dsetID;
  obj_err_t err;

#if defined (EXP_HDF5)
  dsetID = H5Dopen(fileID, dsetStr);
  assert(dsetID >= 0);
#elif defined (EXP_PNETCDF)
  err = ncmpi_inq_varid(fileID, dsetStr, &dsetID);
  assert(err == NC_NOERR);
#endif

  return dsetID;
}


void Release_Dataset_Handle(const obj_handle_t dsetID)
{
  obj_err_t err;

#if defined (EXP_HDF5)
  err = H5Dclose(dsetID);
  assert(err >= 0);
#endif
}


void Query_Dataset_extent(const obj_handle_t fileID,
			  const obj_handle_t dsetID,
			  const int fileFmt,
			  obj_size_t dims[],
			  int *pndims)
{
  obj_handle_t dsetSpace;
  obj_err_t err;
  int dimIDs[5], ndims, hdims, i;

  ndims = (fileFmt == 9 ? 4 : 5);
#if defined (EXP_HDF5)
  dsetSpace = H5Dget_space(dsetID);
  assert(dsetSpace >= 0);
  hdims = H5Sget_simple_extent_dims(dsetSpace, dims, NULL);
  assert(ndims == hdims);
  err = H5Sclose(dsetSpace);
  assert(err >= 0);
#elif defined (EXP_PNETCDF)
  err = ncmpi_inq_vardimid(fileID, dsetID, dimIDs);
  assert(err == NC_NOERR);
  for (i=0; i<ndims; ++i) {
    err = ncmpi_inq_dimlen(fileID, dimIDs[i], &dims[i]);
    assert(err == NC_NOERR);
  }
#endif
  *pndims = ndims;
}


void Read_Dataset(obj_handle_t fileID,
		  obj_handle_t dsetID,
		  obj_size_t dims[],
		  size_t bufsize,
		  double *buf)
{
  obj_err_t err;
  obj_size_t start[5] = {0,0,0,0,0};

#if defined (EXP_HDF5)
  err = H5Dread(dsetID, H5T_NATIVE_DOUBLE, H5S_ALL, 
		H5S_ALL, H5P_DEFAULT, buf);
  assert(err >= 0);
#elif defined (EXP_PNETCDF)
  err = ncmpi_get_vara_all(fileID, dsetID, start, dims, buf, 
                           bufsize, MPI_DOUBLE_PRECISION);
  assert(err == NC_NOERR);
#endif
}


void Read_Attribute(obj_handle_t fileID,
		    obj_handle_t dsetID,
		    const char attName[],
		    double att[])
{
  obj_err_t err;
  obj_handle_t attID;

#if defined (EXP_HDF5)

#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8)
  /* This function is now deprecated by HDF5 */
  attID = H5Aopen_name(dsetID, attName);
#else
  attID = H5Aopen(dsetID, attName, H5P_DEFAULT);
#endif
  assert(attID >= 0);
  err = H5Aread(attID, H5T_NATIVE_DOUBLE, att);
  assert(err >= 0);
  err = H5Aclose(attID);
  assert(err >= 0);

#elif defined (EXP_PNETCDF)
  err = ncmpi_get_att_double(fileID, dsetID, attName, att);
  assert(err == NC_NOERR);
#endif
}


int check_mesh_variables(const double * const buf,
			 const int dsetIndex,
			 const int lenVariables)
{
  int i, fail, val;
  i = 0;
  fail = SIM_SUCCESS;

  while (fail == 0 && i < lenVariables) {

    val = i + dsetIndex + 1; /* expected value */
#ifdef DEBUG_SIM
    printf(" [%s]: (Dataset-%d): actual value %d, expected value %d\n",
	   __FILE__, dsetIndex, (int)buf[i], val);
#endif
    if ((int)buf[i] != val) fail = SIM_FAIL;
    ++i;
  }

  return fail;
}
