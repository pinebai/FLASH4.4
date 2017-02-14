#ifndef IO_H5_XFER_WRAPPER_H
#define IO_H5_XFER_WRAPPER_H

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "io_h5_type.h"
#include "io_h5_xfer.h"
#include "io_h5_type_matched_xfer.h"
#include "io_h5_report_xfer_method.h"
#include "hdf5_flash.h"

int io_h5_xfer_wrapper(const int myPE, const int fileID, const int xferType,
		       const int typeMatchedXfer,
		       const char datasetName[], const int memType,
		       const int memSize[], const int memStart[],
		       const int memCount[], const int diskStart[],
		       const int diskCount[], const int dims, void * pData);
#endif
