#ifndef SIM_EXPECTED_MESH_VALUES_H
#define SIM_EXPECTED_MESH_VALUES_H

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "constants.h"
#include "Flash.h"
#include "io_flash.h"
#include "mangle_names.h"

#if defined (FLASH_IO_EXPERIMENTAL) && defined(FLASH_IO_HDF5)
#define EXP_HDF5
#elif defined (FLASH_IO_EXPERIMENTAL) && defined(FLASH_IO_PNETCDF)
#define EXP_PNETCDF
#endif

#if defined (EXP_HDF5)
#include <hdf5.h>
#include "io_h5_read_file_format.h"
typedef hid_t obj_handle_t;
typedef herr_t obj_err_t;
typedef hsize_t obj_size_t;
#elif defined (EXP_PNETCDF)
#include <pnetcdf.h>
#include "io_ncmpi_read_file_format.h"
typedef int obj_handle_t;
typedef int obj_err_t;
typedef MPI_Offset obj_size_t;
#else
typedef int obj_handle_t;
typedef int obj_err_t;
typedef size_t obj_size_t;
#endif

#define LABEL_LEN 4
#define NUM_ATTR 2
#define SIM_FAIL 1
#define SIM_SUCCESS 0

int FTOC(sim_expected_mesh_values)(const int * const pMyPE,
				   const char * const pLabels);

obj_handle_t Obtain_File_Handle(const char FileName[]);

void Release_File_Handle(obj_handle_t fileID);

int Query_File_Fmt(const int myPE,
		   const obj_handle_t fileID);

obj_handle_t Obtain_Dataset_Handle(const obj_handle_t fileID,
				   const char dsetStr[]);

void Release_Dataset_Handle(const obj_handle_t dsetID);

void Query_Dataset_extent(const obj_handle_t fileID,
			  const obj_handle_t dsetID,
			  const int fileFmt,
			  obj_size_t dims[],
			  int *pndims);

void Read_Dataset(obj_handle_t fileID,
		  obj_handle_t dsetID,
		  obj_size_t dims[],
		  size_t bufsize,
		  double *buf);

void Read_Attribute(obj_handle_t fileID,
		    obj_handle_t dsetID,
		    const char attName[],
		    double att[]);

int check_mesh_variables(const double * const buf,
			 const int dsetIndex,
			 const int lenVariables);
#endif
