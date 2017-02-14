#ifndef IO_NCMPI_XFER_MESH_DATASET_H
#define IO_NCMPI_XFER_MESH_DATASET_H

#include <pnetcdf.h>
#include <assert.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "mangle_names.h"
#include "io_use_grid_mpi_types.h"
#include "io_ncmpi_nonblocking.h"
#include "io_ncmpi_name_to_id.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

/* Test for Pnetcdf >= 1.2.0 */
#if ((PNETCDF_VERSION_MAJOR > 1) ||					\
     (PNETCDF_VERSION_MAJOR == 1 && PNETCDF_VERSION_MINOR >= 2))

#define flash_iput_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iput_vara(a,b,c,d,e,f,g,h)
#define flash_iget_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iget_vara(a,b,c,d,e,f,g,h)

#else

/* We can detect old versions of pnetcdf as there is no
   PNETCDF_VERSION_MAJOR version macro in pnetcdf.h */
#ifdef PNETCDF_VERSION_MAJOR
/* Version 1.1.1 */
#define flash_iput_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iput_vara_all(a,b,c,d,e,f,g,h)
#define flash_iget_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iget_vara_all(a,b,c,d,e,f,g,h)
#else
/* Version <= 1.1.0 */
#define flash_iput_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iput_vara(a,b,c,d,e,f,g,h)
#define flash_iget_vara(a,b,c,d,e,f,g,h) \
  ncmpi_iget_vara(a,b,c,d,e,f,g,h)
#endif

#endif

void io_ncmpi_xfer_mesh_dataset(const int myPE,
				const int fileID,
				const int xferType,
				const int fileType,
				const int fileFmt,
				const int gridDataStruct,
				const int numDataDims,
				const int nonBlocking,
				const int globalOffset[],
				const int localSubSize[],
				const char datasetName[],
				double * pData);
#endif
