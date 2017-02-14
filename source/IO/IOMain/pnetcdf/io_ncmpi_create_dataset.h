#ifndef IO_NCMPI_CREATE_DATASET_H
#define IO_NCMPI_CREATE_DATASET_H

#include "io_ncmpi_type.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <pnetcdf.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

void io_ncmpi_create_dataset(const int myPE,
			     const int fileID,
			     const int diskType,
			     const int dims,
			     const int dimIDs[],
			     const char datasetName[]);

#endif
