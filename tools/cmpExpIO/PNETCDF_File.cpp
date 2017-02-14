#include "PNETCDF_File.hpp"

//Constructor reads the list of variables from the file.
PNETCDF_File::PNETCDF_File(const WorkCoordinator *workCoordinator, 
			   const std::string &fileName,
			   const std::string &meshName,
			   const bool &useCollectiveIO,
			   const bool &standardFlashFile)
{
  int err, i, varID;
  const char * pQueryVar = NULL;

  m_useCollectiveIO = useCollectiveIO;
  m_standardFlashFile = standardFlashFile;
  m_meshName = meshName;

  err = ncmpi_open(MPI_COMM_WORLD, fileName.c_str(), NC_NOWRITE, 
		   MPI_INFO_NULL, &m_ncID);
  assert(err == NC_NOERR);

  if (!m_useCollectiveIO) {
    err = ncmpi_begin_indep_data(m_ncID);
    assert(err == NC_NOERR);
  }

  m_variableMap = ReadMeshVarNamesFromFile(m_ncID, "unknown");
  if (m_variableMap.size() == 0) {
    std::cout << "Exiting because there are zero mesh variables in file " 
	      << fileName << std::endl;
    exit(EXIT_SUCCESS);
  }

  /* Find the extent of each dimension of the mesh variable(s) */
  pQueryVar = m_standardFlashFile ?
    m_variableMap.begin()->first.c_str() : m_meshName.c_str();

  err = ncmpi_inq_varid(m_ncID, pQueryVar, &varID);
  assert(err == NC_NOERR);

  err = ncmpi_inq_varndims(m_ncID, varID, &m_numDims);
  assert(err == NC_NOERR && (m_numDims == 4 || m_numDims == 5));

  err = ncmpi_inq_vardimid(m_ncID, varID, m_dimIDs);
  assert(err == NC_NOERR);

  for (i=0; i<m_numDims; ++i) {
    err = ncmpi_inq_dimlen(m_ncID, m_dimIDs[i], &m_dims[i]);
    assert(err == NC_NOERR);
  }

  /* Figure out which processor takes care of what blocks (m_dims[0]=blocks) */
  workCoordinator->GetWorkPortion(m_dims[0], &m_myBlockStart, &m_myBlockCount);
  assert(m_myBlockCount <= m_dims[0]);
  assert(m_myBlockStart >= 0 && m_myBlockStart < m_dims[0]);

  m_numGridPoints = m_myBlockCount * m_dims[1] * m_dims[2] * m_dims[3];

  if (m_myBlockCount > 0) {
    m_data = new double[m_numGridPoints];
    assert(NULL != m_data);
  }
}


PNETCDF_File::~PNETCDF_File()
{
  int err;

  if (!m_useCollectiveIO) {
    err = ncmpi_end_indep_data(m_ncID);
    assert(err == NC_NOERR);
  }

  err = ncmpi_close(m_ncID);
  assert(err == NC_NOERR);

  delete[] m_data;
}


std::vector<std::string> PNETCDF_File::GetAllVariableNames() const
{
  std::vector<std::string> v;
  for(MapType::const_iterator it = m_variableMap.begin(); 
      it != m_variableMap.end(); ++it) {
    v.push_back(it->first);
  }
  return v;
}


size_t PNETCDF_File::GetNumberDataElements() const
{
  return m_numGridPoints;
}


const double * PNETCDF_File::GetVariableFromFile(const std::string &varName) const
{
  MapType::const_iterator it = m_variableMap.find(varName);
  MPI_Offset start[MAXDIMS], stride[MAXDIMS], count[MAXDIMS];
  int err, varID;
  const char *pVarToRead = NULL;

  /* Return a null pointer if the file does not contain the variable name */
  if (it == m_variableMap.end()) return NULL;
  pVarToRead = m_standardFlashFile ? varName.c_str() : m_meshName.c_str();

  err = ncmpi_inq_varid(m_ncID, pVarToRead, &varID);
  assert(err == NC_NOERR);

  /* Each variable is given its own space in the file.  So we just need to 
     distribute the data across processors. */
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


  if (m_useCollectiveIO) {
    /* If we are not in the default collective mode when we call
       ncmpi_get_vars_all then we get an error return code of NC_EINDEP(-27) */
    err = ncmpi_get_vars_all(m_ncID, varID, start, count, stride, &m_data[0], 
			     m_numGridPoints, MPI_DOUBLE_PRECISION);
  } else {
    /* If we are not in independent mode when we call
       ncmpi_get_vars then we get an error return code of NC_ENOTINDEP(-26) */
    err = ncmpi_get_vars(m_ncID, varID, start, count, stride, &m_data[0], 
			 m_numGridPoints, MPI_DOUBLE_PRECISION);
  }
  assert(err == NC_NOERR);

  return &m_data[0];
}


MapType PNETCDF_File::ReadMeshVarNamesFromFile(const int &dataFileID,
					       const char meshPrefix[])
{
  MapType stringMap;
  MPI_Offset attLen;
  int err, i, nVar, lenStrID;
  char attNameNumVar[BIG_STRING_LEN], attNameStr[BIG_STRING_LEN], 
    buf[BIG_STRING_LEN], strID[BIG_STRING_LEN];


  /* Construct the name of the attribute containing the number
     of mesh variables */
  strncpy(attNameNumVar, meshPrefix, BIG_STRING_LEN);
  strncat(attNameNumVar, "_names", BIG_STRING_LEN);
  err = ncmpi_get_att_int(dataFileID, NC_GLOBAL, attNameNumVar, &nVar);
  assert(err == NC_NOERR && nVar > 0);


  strncpy(attNameStr, meshPrefix, BIG_STRING_LEN);
  strncat(attNameStr, "_", BIG_STRING_LEN);
  lenStrID = strlen(attNameStr);
  for (i=0; i<nVar; ++i) {

    /* Construct the pnetcdf attribute name */
    snprintf(strID, BIG_STRING_LEN, "%d", i);
    strncpy(&attNameStr[lenStrID], strID, BIG_STRING_LEN);

    /* Figure out attribute length and ensure it is < BIG_STRING_LEN */
    err = ncmpi_inq_attlen(dataFileID, NC_GLOBAL, attNameStr, &attLen);
    assert(err == NC_NOERR && attLen < BIG_STRING_LEN);
    
    /* Read attribute from file */
    err = ncmpi_get_att_text(dataFileID, NC_GLOBAL, attNameStr, buf);
    assert(err == NC_NOERR);
    buf[attLen] = '\0';
    
    /* Store in a map of FLASH variable names */
    stringMap[buf] = i;
  }
  return stringMap;
}
