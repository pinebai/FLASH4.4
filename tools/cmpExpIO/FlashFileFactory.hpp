#ifndef FLASHFILEFACTORY_H
#define FLASHFILEFACTORY_H

#include "mpi.h"
#include "FlashFile.hpp"
#include "WorkCoordinator.hpp"

#ifndef NO_HDF5
#include "HDF5_File.hpp"
#endif

#ifndef NO_NCDF
#include "PNETCDF_File.hpp"
#endif


class FlashFileFactory {
private:
  WorkCoordinator *m_workCoordinator;
  int m_myPE, m_numProcs;
  bool m_useCollectiveIO;

public:
  FlashFileFactory(const int &myPE,
		   const int &numProcs,
		   const bool &useCollectiveIO);
  ~FlashFileFactory();

  FlashFile * GetFlashFileInstance(const std::string &fileName) const;

#ifndef NO_HDF5
  bool IsFileHDF5(const std::string &fileName) const;
  int GetHDF5FileFormatVersion(const std::string &fileName) const;
#endif

#ifndef NO_NCDF
  bool IsFilePNETCDF(const std::string &fileName) const;
  int GetPNETCDFFileFormatVersion(const std::string &fileName) const;
#endif
};

#endif
