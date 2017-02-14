#include "HDF5_File.hpp"

//Constructor reads the list of variables from the file.
HDF5_File::HDF5_File(const WorkCoordinator *workCoordinator, 
		     const std::string &fileName,
		     const std::string &meshName,
		     const bool &useCollectiveIO,
		     const bool &standardFlashFile)
{
  hid_t singleVarDset, singleVarDataSpace;
  herr_t err;
  int ndims;
  hsize_t start[4], stride[4], count[4];  
  size_t numVariables;
  const char * pQueryVar = NULL;
  H5FD_mpio_xfer_t xfer_mode;

  m_useCollectiveIO = useCollectiveIO;
  m_standardFlashFile = standardFlashFile;
  m_meshName = meshName;
  m_data = NULL;

  /* Set up file access property list with parallel I/O access */
  hid_t fileAccessListID = H5Pcreate(H5P_FILE_ACCESS);
  assert(fileAccessListID != -1);

  err = H5Pset_fapl_mpio(fileAccessListID, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(err >= 0);

  /* Open file collectively and release property list identifier. */
  m_file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, fileAccessListID);
  assert(m_file >= 0);
  err = H5Pclose(fileAccessListID);
  assert(err >= 0);

  /* Create the property list for collective or independent transfers */
  m_h5_xferList = H5Pcreate(H5P_DATASET_XFER);
  assert(m_h5_xferList != -1);

  xfer_mode = useCollectiveIO ? H5FD_MPIO_COLLECTIVE : H5FD_MPIO_INDEPENDENT;
  err = H5Pset_dxpl_mpio(m_h5_xferList, xfer_mode);
  assert(err >= 0);

  /* Read the name of each FLASH variable in the mesh dataset */
  m_variableMap = ReadStringsFromFile(m_file, "unknown names", m_h5_xferList);
  numVariables = m_variableMap.size(); /* Must be at least 1 variable */
  if (numVariables == 0) {
    std::cout << "Exiting because there are zero mesh variables in file " 
	      << fileName << std::endl;
    exit(EXIT_SUCCESS);
  }

  /* Find the extent of each dimension of the mesh variable(s) */
  pQueryVar = m_standardFlashFile ?
    m_variableMap.begin()->first.c_str() : m_meshName.c_str();

  singleVarDset = H5Dopen(m_file, pQueryVar);
  assert(singleVarDset >= 0);

  singleVarDataSpace = H5Dget_space(singleVarDset);
  assert(singleVarDataSpace >= 0);

  ndims = H5Sget_simple_extent_dims(singleVarDataSpace, m_dims, NULL);
  assert(ndims == 4 || ndims == 5);

  err = H5Sclose(singleVarDataSpace);
  assert(err >= 0);

  err = H5Dclose(singleVarDset);
  assert(err >= 0);

  /* Figure out which processor takes care of what blocks (m_dims[0]=blocks) */
  workCoordinator->GetWorkPortion(m_dims[0], &m_myBlockStart, &m_myBlockCount);
  assert(m_myBlockCount <= m_dims[0]);
  assert(m_myBlockStart >= 0 && m_myBlockStart < m_dims[0]);

  m_numGridPoints = m_myBlockCount * m_dims[1] * m_dims[2] * m_dims[3];

  if (m_myBlockCount > 0) {
    m_data = new double[m_numGridPoints];
    assert(NULL != m_data);

    /* Define the memory dataspace. */
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = m_myBlockCount;
    count[1] = m_dims[1];
    count[2] = m_dims[2];
    count[3] = m_dims[3];

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;
    stride[3] = 1;

    m_memSpace = H5Screate_simple(4, count, NULL);
    assert(m_memSpace >= 0);
    err = H5Sselect_hyperslab(m_memSpace, H5S_SELECT_SET, 
			      start, stride, count, NULL);
    assert(err >= 0);
  } else {

    count[0] = 1;
    count[1] = 1;
    count[2] = 1;
    count[3] = 1;

    m_memSpace = H5Screate_simple(4, count, NULL);
    assert(m_memSpace >= 0);

    err = H5Sselect_none(m_memSpace);
    assert(err >= 0);
  }
}


HDF5_File::~HDF5_File()
{
  herr_t err;

  delete [] m_data;

  err = H5Pclose(m_h5_xferList);
  assert(err >= 0);

  err = H5Sclose(m_memSpace);
  assert(err >= 0);

  err = H5Fclose(m_file);
  assert(err >= 0);
}


std::vector<std::string> HDF5_File::GetAllVariableNames() const
{
  std::vector<std::string> v;
  for(MapType::const_iterator it = m_variableMap.begin(); 
      it != m_variableMap.end(); ++it) {
    v.push_back(it->first);
  }
  return v;
}


size_t HDF5_File::GetNumberDataElements() const
{
  return m_numGridPoints;
}


const double * HDF5_File::GetVariableFromFile(const std::string &varName) const
{
  MapType::const_iterator it = m_variableMap.find(varName);
  hsize_t start[MAXDIMS], stride[MAXDIMS], count[MAXDIMS];
  hid_t datasetID, dataSpace;
  herr_t err;
  const char *pVarToRead = NULL;

  /* Return a null pointer if the file does not contain the variable name */
  if (it == m_variableMap.end()) return NULL;
  pVarToRead = m_standardFlashFile ? varName.c_str() : m_meshName.c_str();

  /* Open single variable dataset (only appears in regular checkpoint file). */
  datasetID = H5Dopen(m_file, pVarToRead);
  assert(datasetID >= 0);
  dataSpace = H5Dget_space(datasetID);
  assert(dataSpace >= 0);

  /* Define a hyperslab that selects a single variable from the 
     data on disk */
  if (m_myBlockCount > 0) {
    start[0] = m_myBlockStart;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = m_myBlockCount;
    count[1] = m_dims[1];
    count[2] = m_dims[2];
    count[3] = m_dims[3];

    stride[0] = 1;
    stride[1] = 1;
    stride[2] = 1;
    stride[3] = 1;

    if (!m_standardFlashFile) {
      start[4] = it->second;
      count[4] = 1;
      stride[4] = m_dims[4];
    }

    err = H5Sselect_hyperslab(dataSpace, H5S_SELECT_SET,
			      start, stride, count, NULL);
    assert(err >= 0);
  } else {
    err = H5Sselect_none(dataSpace);
    assert(err >= 0);
  }

  err = H5Dread(datasetID, H5T_NATIVE_DOUBLE, m_memSpace, 
                dataSpace, m_h5_xferList, &m_data[0]);
  assert(err >= 0);

  /* Close the dataspace. */
  err = H5Sclose(dataSpace);
  assert(err >= 0);

  /* Close the dataset. */
  err = H5Dclose(datasetID);
  assert(err >= 0);

  return &m_data[0];
}


/* See HDF5 sample application: h5ex_t_string.c */
MapType HDF5_File::ReadStringsFromFile(const hid_t &dataFile,
				       const char dataSetName[],
				       const hid_t &xferList)
{
  MapType stringMap;
  hsize_t dims[2];
  size_t strSize, i;
  hid_t dataset_id, dataspace_id, datatype_id, memDatatype;
  herr_t err;
  int ndims;
  char **rdata;     


  dataset_id = H5Dopen(dataFile, dataSetName);
  assert(dataset_id >= 0);

  dataspace_id = H5Dget_space(dataset_id);
  assert(dataspace_id >= 0);


  /* ----------------------------------------------------------------
     Query the string datatype stored in the file and then create 
     an identical datatype in memory.
     ---------------------------------------------------------------- */
  datatype_id = H5Dget_type(dataset_id);
  assert(datatype_id >= 0);

  strSize = H5Tget_size(datatype_id);
  assert(strSize != 0);
  strSize++; /* Space for null terminator */

  memDatatype = H5Tcopy(H5T_C_S1);
  assert(memDatatype >=0);

  err = H5Tset_size(memDatatype, strSize);
  assert(err >= 0);

  ndims = H5Sget_simple_extent_ndims(dataspace_id);
  assert(ndims == 2); /* We expect array of strings in the file */

  ndims = H5Sget_simple_extent_dims (dataspace_id, dims, NULL);
  assert(ndims == 2);


  /* ----------------------------------------------------------------
     Allocate enough memory so we can do a single read from the file.
     Copy this data into an STL map of strings and then free memory.
     ---------------------------------------------------------------- */
  rdata = (char **) malloc(dims[0] * sizeof (char *));
  assert(rdata != NULL);
  rdata[0] = (char *) malloc(dims[0] * strSize * sizeof (char));
  assert(rdata[0] != NULL);
  for (i=1; i<dims[0]; i++)
    rdata[i] = rdata[0] + i * strSize;

  err = H5Dread(dataset_id, memDatatype, H5S_ALL, H5S_ALL, xferList, rdata[0]);
  assert(err >= 0);

  for (i=0; i<dims[0]; i++) {
    stringMap[rdata[i]] = i;
  }

  free(rdata[0]);
  free(rdata);


  /* ----------------------------------------------------------------
     Free all HDF5 resources
     ---------------------------------------------------------------- */
  err = H5Tclose(memDatatype);
  assert(err >= 0);

  err = H5Sclose(dataspace_id);
  assert(err >= 0);

  err = H5Dclose(dataset_id);
  assert(err >= 0);


  return stringMap;
}
