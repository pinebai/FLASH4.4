#ifndef IO_H5_XFER_MESH_DATASET_H
#define IO_H5_XFER_MESH_DATASET_H

#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "hdf5_flash.h"
#include <hdf5.h>
#include <assert.h>
#include <string.h>
#include "io_h5_report_xfer_method.h"

void io_h5_xfer_mesh_dataset(const int myPE,
			     const int fileID,
			     const int xferType,
			     const int fileFmt,
			     const int fileType,
			     const int gridDataStruct,
			     const int numFileDims,
			     const int numGridVar,
			     const int outputVarOffset[],
			     const int globalOffset[],
			     const int localSize[],
			     const int localSubSize[],
			     const int localOffset[],
			     const char datasetName[],
			     double * pData);
#endif
