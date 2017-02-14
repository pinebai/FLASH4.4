#ifndef IO_CREATE_DATASET_H
#define IO_CREATE_DATASET_H

#include "mangle_names.h"
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

#if defined(FLASH_IO_PNETCDF)
#include "io_ncmpi_create_dataset.h"
#endif
#if defined(FLASH_IO_HDF5)
#include "io_h5create_dataset.h"
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

void FTOC(io_create_dataset)(const int * const pMyPE,
			     const int * const pFileID,
			     const int * const pLibType,
			     const int * const pDiskType,
			     const int * const pDims,
			     const int dimIDs[],
			     const char datasetName[],
			     const int * const pDsetNameLen);
#endif
