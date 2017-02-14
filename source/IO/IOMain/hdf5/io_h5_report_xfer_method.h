#ifndef IO_H5_REPORT_XFER_METHOD_H
#define IO_H5_REPORT_XFER_METHOD_H

#include <hdf5.h>
#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"

/* Useful macro which is only defined in HDF5 >= 1.8.7 (see H5public.h) */
#ifndef H5_VERSION_GE
#define H5_VERSION_GE(Maj,Min,Rel) \
  (((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR==Min) && (H5_VERS_RELEASE>=Rel)) || \
  ((H5_VERS_MAJOR==Maj) && (H5_VERS_MINOR>Min)) || \
   (H5_VERS_MAJOR>Maj))
#endif

#define BUFSIZE 1000

int io_h5_report_xfer_method(const int myPE, const hid_t hXferList,
			     const char datasetName[]);

#endif
