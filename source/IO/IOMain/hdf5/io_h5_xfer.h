#ifndef IO_H5_XFER_H
#define IO_H5_XFER_H

#include <hdf5.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

int io_h5_xfer(const int myPE, const hid_t fileID, const int xferType,
	       const hid_t hXferList, const char datasetName[],
	       const hid_t hMemType, const hsize_t hMemSize[],
	       const hsize_t hMemStart[], const hsize_t hMemCount[],
	       const hsize_t hDiskStart[], const hsize_t hDiskCount[],
	       const int dims, void * pData);
#endif
