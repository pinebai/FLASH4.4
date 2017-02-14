#ifndef IO_H5CREATE_DATASET_H
#define IO_H5CREATE_DATASET_H

#include <hdf5.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "io_h5_type.h"

void io_h5create_dataset(const int myPE,
			 const int fileID,
			 const int diskType,
			 const int dims,
			 const int diskSize[],
			 const char dataSetName[]);

#endif
