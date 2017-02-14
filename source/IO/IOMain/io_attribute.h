#ifndef IO_ATTRIBUTE_H
#define IO_ATTRIBUTE_H

#include "mangle_names.h"
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

#if defined(FLASH_IO_PNETCDF)
#include "io_ncmpi_attribute.h"
#endif
#if defined(FLASH_IO_HDF5)
#include "io_h5_attribute.h"
#endif

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

void FTOC(io_attribute_create)(const int * const pMyPE,
			       const int * const pFileID,
			       const int * const pLibType,
			       const int * const pDiskType,
			       const int * const pDims,
			       const int datasetSize[],
			       const char datasetName[],
			       const int * const pDsetNameLen,
			       const char attDatasetName[],
  			       const int * const pAttNameLen);

void FTOC(io_attribute_write)(const int * const pMyPE,
			      const int * const pFileID,
			      const int * const pLibType,
			      const int * const pMemType,
			      const char datasetName[],
			      const int * const pDsetNameLen,
			      const char attDatasetName[],
			      const int * const pAttNameLen,
			      const void * const pData);
#endif
