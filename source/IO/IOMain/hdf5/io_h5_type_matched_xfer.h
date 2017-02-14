#ifndef IO_H5_TYPE_MATCHED_XFER_H
#define IO_H5_TYPE_MATCHED_XFER_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "io_h5_type.h"
#include "io_h5_xfer.h"
#include "io_repack_data.h"
#include "hdf5_flash.h"

int
io_h5_type_matched_xfer(const int myPE, const hid_t hFileID, const int xferType,
			const hid_t hXferList, const char datasetName[],
			const hid_t hMemType, const hsize_t hMemSize[],
			const hsize_t hMemStart[], const hsize_t hMemCount[],
			const hsize_t hDiskStart[], const hsize_t hDiskCount[],
			const int dims, void * pData);

int
io_h5_use_extra_float_buffer(const int myPE, const hid_t hFileID,
			     const char datasetName[], const hid_t hMemType);

#endif
