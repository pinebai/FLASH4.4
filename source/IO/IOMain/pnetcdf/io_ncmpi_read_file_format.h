#ifndef IO_NCMPI_READ_FILE_FORMAT_H
#define IO_NCMPI_READ_FILE_FORMAT_H

#include <mpi.h>
#include <pnetcdf.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "mangle_names.h"

#ifdef USE_IO_C_INTERFACE
#ifdef FTOC
#undef FTOC
#endif
#define FTOC(x) x
#endif

void FTOC(io_ncmpi_read_file_format)(const int * const pMyPE,
				     const int * const pFileID,
				     int *pFileFormat);

#endif
