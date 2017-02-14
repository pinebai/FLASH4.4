#ifndef PNETCDF_FILE_H
#define PNETCDF_FILE_H

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <map>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "pnetcdf.h"
#include "flash_types.hpp"
#include "FlashFile.hpp"
#include "WorkCoordinator.hpp"

class PNETCDF_File : public FlashFile {

private:
  MapType m_variableMap;
  mutable double *m_data;
  size_t m_numGridPoints;
  MPI_Offset m_dims[MAXDIMS];
  int m_dimIDs[MAXDIMS], m_ncID, m_myBlockStart, m_myBlockCount, m_numDims;
  bool m_standardFlashFile, m_useCollectiveIO;
  std::string m_meshName;
  MapType ReadMeshVarNamesFromFile(const int &dataFileID,
				   const char meshPrefix[]);

public:
  PNETCDF_File(const WorkCoordinator *workCoordinator, 
	       const std::string &fileName,
	       const std::string &meshName,
	       const bool &useCollectiveIO,
	       const bool &standardFlashFile);
  ~PNETCDF_File();

  std::vector<std::string> GetAllVariableNames() const;
  size_t GetNumberDataElements() const;
  const double * GetVariableFromFile(const std::string &varName) const;
};

#endif
