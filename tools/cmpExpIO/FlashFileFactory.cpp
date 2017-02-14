#include "FlashFileFactory.hpp"

FlashFileFactory::FlashFileFactory(const int &myPE,
				   const int &numProcs,
				   const bool &useCollectiveIO)
{
  m_myPE = myPE;
  m_numProcs = numProcs;
  m_useCollectiveIO = useCollectiveIO;
  m_workCoordinator = new WorkCoordinator(myPE, numProcs);
}

FlashFileFactory::~FlashFileFactory()
{
  delete m_workCoordinator;
}

FlashFile * FlashFileFactory::GetFlashFileInstance(
  const std::string &fileName) const
{
  int versionID;
  std::string errMsg, infoMsg;
  bool standardFlashFile;

  if (m_myPE == 0) {
    std::cout << "Testing type of file " << fileName << std::endl;
  }

#ifndef NO_HDF5
  if (IsFileHDF5(fileName)) {
    
    infoMsg = "File " + fileName + " is hdf5 format";
    versionID = GetHDF5FileFormatVersion(fileName);

    if (versionID == 9) {
      infoMsg += " and has regular FLASH hdf5 format";
      standardFlashFile = true;
    } else if (versionID == 10) {
      infoMsg += " and has new FLASH hdf5 format";
      standardFlashFile = false;
    } else {
      errMsg = "File " + fileName + " has unrecognised hdf5 version number";
      std::cerr << errMsg << std::endl;
      exit(EXIT_FAILURE);
    }

    if (m_myPE == 0) {
      std::cout << infoMsg << std::endl;
    }
    return new HDF5_File(m_workCoordinator, fileName, "unknown",
			 m_useCollectiveIO, standardFlashFile);
  }
#endif

#ifndef NO_NCDF
  if (IsFilePNETCDF(fileName)) {

    infoMsg = "File " + fileName + " is pnetcdf format";
    versionID = GetPNETCDFFileFormatVersion(fileName);

    if (versionID == 9) {
      standardFlashFile = true;
      infoMsg += " and has regular FLASH pnetcdf format";
    } else if (versionID == 10) {
      standardFlashFile = false;
      infoMsg += " and has new FLASH pnetcdf format";
    } else {
      errMsg = "File " + fileName + " has unrecognised pnetcdf version number";
      std::cerr << errMsg << std::endl;
      exit(EXIT_FAILURE);
    }

    if (m_myPE == 0) {
      std::cout << infoMsg << std::endl;
    }
    return new PNETCDF_File(m_workCoordinator, fileName, "unknown",
			    m_useCollectiveIO, standardFlashFile);
  }
#endif

  std::cerr << "File " << fileName << " not recognised\n";
  exit(EXIT_FAILURE);
  return NULL;
}


#ifndef NO_HDF5
bool FlashFileFactory::IsFileHDF5(const std::string &fileName) const
{
  /* Special case for when we use bglockless on BG/P.  ROMIO 
     understands the prefix "bglockless:" but H5Fis_hdf5 does not */
  htri_t status;
  std::string pureFileName;
  const std::string bglocklessStr = "bglockless:";
  std::string::size_type loc = fileName.find(bglocklessStr, 0);
  if(loc == 0) {
    /* Remove substring only when it begins at the first character */
    pureFileName = fileName.substr(bglocklessStr.length());
  } else {
    pureFileName = fileName;
  }

  status = H5Fis_hdf5(pureFileName.c_str());
  assert(status >= 0);
  return (status > 0); /* A positive value for status indicates a HDF5 file */
}

int FlashFileFactory::GetHDF5FileFormatVersion(
  const std::string &fileName) const
{
  int versionNum;
  herr_t err;
  hid_t file, fileAccessListID, dset, dtypeInt;

  /* Set up file access property list with parallel I/O access */
  fileAccessListID = H5Pcreate(H5P_FILE_ACCESS);
  assert(fileAccessListID != -1);

  err = H5Pset_fapl_mpio(fileAccessListID, MPI_COMM_WORLD, MPI_INFO_NULL);
  assert(err >= 0);

  /* Open file collectively and release property list identifier. */
  file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, fileAccessListID);
  assert(file >= 0);
  err = H5Pclose(fileAccessListID);
  assert(err >= 0);


  dset = H5Dopen(file, "sim info");
  assert(dset >= 0);
  dtypeInt = H5Tcreate(H5T_COMPOUND, sizeof(int));    
  assert(dtypeInt >= 0);
  err = H5Tinsert(dtypeInt, "file format version", 0, H5T_NATIVE_INT);
  assert(err >= 0);


  err = H5Dread(dset, dtypeInt, H5S_ALL, H5S_ALL, H5P_DEFAULT, &versionNum);
  assert(err >= 0);


  err = H5Tclose(dtypeInt);
  assert(err >= 0);
  err = H5Dclose(dset);
  assert(err >= 0);
  err = H5Fclose(file);
  assert(err >= 0);

  return versionNum;
}
#endif


#ifndef NO_NCDF
bool FlashFileFactory::IsFilePNETCDF(const std::string &fileName) const
{
  int err, format;
  /* I ASSUME THAT PNETCDF DOES NOT MODIFY fileName */
  char *pFileName;
  pFileName = const_cast<char*>(fileName.c_str());
  err = ncmpi_inq_file_format(pFileName, &format);
  assert (err == NC_NOERR || err == NC_ENOTNC);
  return (format != NC_FORMAT_UNKNOWN);
}

int FlashFileFactory::GetPNETCDFFileFormatVersion(
  const std::string &fileName) const
{
  int versionNum, err, ncID;
  err = ncmpi_open(MPI_COMM_WORLD, fileName.c_str(), NC_NOWRITE, 
		   MPI_INFO_NULL, &ncID);

  assert(err == NC_NOERR);

  err = ncmpi_get_att_int(ncID, NC_GLOBAL, "file_format_version", &versionNum);
  assert(err == NC_NOERR);

  err = ncmpi_close(ncID);
  assert(err == NC_NOERR);

  return versionNum;
}
#endif
