#ifndef IO_NCMPI_NAME_TO_ID_H
#define IO_NCMPI_NAME_TO_ID_H

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

void FTOC(io_ncmpi_name_to_id)(const int * const pFileID,
			       const char name[],
			       const int * const pLen,
			       int *pVarID);

void io_ncmpi_name_to_id_c(const int fileID,
			   const char name[],
			   int *pVarID);

#endif
