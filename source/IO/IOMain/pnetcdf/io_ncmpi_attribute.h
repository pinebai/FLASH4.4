#ifndef IO_NCMPI_ATTRIBUTE_H
#define IO_NCMPI_ATTRIBUTE_H

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "io_ncmpi_type.h"
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

void io_ncmpi_attribute_create(const int myPE,
			       const int fileID,
			       const int diskType,
			       const int dims,
			       const int datasetSize[],
			       const char datasetName[],
			       const char attDatasetName[]);

void io_ncmpi_attribute_write(const int myPE,
			      const int fileID,
			      const int memType,
			      const char datasetName[],
			      const char attDatasetName[],
			      const void * const pData);
#endif
