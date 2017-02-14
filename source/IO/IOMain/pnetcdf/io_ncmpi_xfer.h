#ifndef IO_NCMPI_XFER_H
#define IO_NCMPI_XFER_H

#include <mpi.h>
#include <pnetcdf.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

int io_ncmpi_xfer(const int myPE, const int fileID, const int xferType,
		  const int varID, const MPI_Datatype memType,
		  const MPI_Offset memCountScalar, const MPI_Offset diskStart[],
		  const MPI_Offset diskCount[], const int dims, void * pData);
#endif
