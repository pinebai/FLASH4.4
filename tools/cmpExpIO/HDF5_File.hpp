#ifndef HDF5_FILE_H
#define HDF5_FILE_H

#include <cstdlib>
#include <cassert>
#include <string>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include "hdf5.h"
#include "flash_types.hpp"
#include "FlashFile.hpp"
#include "WorkCoordinator.hpp"

class HDF5_File : public FlashFile {
private:
  MapType m_variableMap;
  mutable double *m_data;
  hid_t m_h5_xferList, m_file, m_memSpace;
  size_t m_numGridPoints;
  hsize_t m_dims[MAXDIMS];
  int m_myBlockStart, m_myBlockCount;
  bool m_standardFlashFile, m_useCollectiveIO;
  std::string m_meshName;

public:
  HDF5_File(const WorkCoordinator *workCoordinator, 
	    const std::string &fileName,
	    const std::string &meshName,
	    const bool &useCollectiveIO,
	    const bool &standardFlashFile);
  ~HDF5_File();

  std::vector<std::string> GetAllVariableNames() const;
  size_t GetNumberDataElements() const;
  const double * GetVariableFromFile(const std::string &varName) const;
  MapType ReadStringsFromFile(const hid_t &dataFile,
			      const char dataSetName[],
			      const hid_t &xferList);
};

#endif
