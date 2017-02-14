#ifndef IO_NCMPI_XFER_WRAPPER_H
#define IO_NCMPI_XFER_WRAPPER_H

#include <mpi.h>
#include <pnetcdf.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "io_mpi_type.h"
#include "io_ncmpi_xfer.h"

int io_ncmpi_xfer_wrapper(const int myPE, const int fileID, const int xferType,
			  const char datasetName[], const int memType,
			  const int memSize[], const int memStart[],
			  const int memCount[], const int diskStart[],
			  const int diskCount[], const int dims, void * pData);
#endif
