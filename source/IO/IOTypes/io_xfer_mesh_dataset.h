#ifndef IO_XFER_MESH_DATASET_H
#define IO_XFER_MESH_DATASET_H

#include "constants.h"
#include "mangle_names.h"
#include "io_flash.h"
#include <assert.h>

#ifdef FLASH_IO_PNETCDF
#include "io_ncmpi_xfer_mesh_dataset.h"
#endif

#ifdef FLASH_IO_HDF5
#include "io_h5_xfer_packed_mesh_dataset.h"
#include "io_h5_xfer_mesh_dataset.h"
#endif

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

void FTOC(io_xfer_mesh_dataset)(const int * const pMyPE,
				const int * const pFileID,
				const int * const pLibType,
				const int * const pXferType,
				const int * const pFileType,
				const int * const pFileFmt,
				const int * const pGridDataStruct,
				const int * const pNumDataDims,
				const int * const pNumGridVars,
				const int * const pNonBlocking,
				const int * const pPrePackData,
				const int diskOffset[],
				const int memSize[],
				const int memSubSize[],
				const int memOffset[],
				const int memVarOffset[],
				const char datasetName[],
				const int * const pDsetNameLen,
				void * pData);

#endif
