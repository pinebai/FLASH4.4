
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/parallel/io_h5read_generic_real_arr_
 *
 * NAME
 *
 *  io_h5read_generic_real_arr_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_generic_real_arr_(file_identifier, generic_arr, local_size, total_size,
 *                        global_offset, dataset_name, name_len)
 *
 *  void io_h5read_generic_real_arr_(hid_t*, int*, int*, int*, int*, char*, int*)
 *
 *
 * DESCRIPTION
 *
 *  Read in a single processors's portion of the data
 *  (generic_arr) from a FLASH HDF5 checkpoint file.  Given the HDF5 file handle 
 *  (file_identifier), the offset to start reading at (global_offset) read in
 *  local_size worth of data.  
 *
 *
 * ARGUMENTS
 *
 *  file_identifier          the HDF5 file handle for the current file (as
 *                           returned by io_h5open_file_for_read_)
 *
 *  generic_arr                  the refinement level data (returned)
 *
 *  local_size             the number of blocks worth of data to read
 *
 *  total_size             the total number of blocks stored in the file
 *                           (as returned from io_h5readHeaderInfo_)
 *
 *  global_offset            the starting position of the current processor's
 *                           data (specified in blocks)
 * 
 *  dataset_name          name of the dataset to read in
 *
 *  name_len              length of the dataset name. makes things much easier
 *                        when trying to convert from c to fortran strings
 *
 ****/

void FTOC(io_h5read_generic_real_arr)(hid_t* file_identifier,
				     int* generic_arr,
				     int* local_size,
				     int* total_size,
				     int* global_offset,
				     char* dataset_name, 
				     int* name_len)
{
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_1d;

  hsize_t start_1d;
  hsize_t stride_1d, count_1d;


  char* dataset_name_new;
  
  dataset_name_new = (char *) malloc((*name_len) + 1 * sizeof(char)); 
  
  /* copy the dataset_name into a c string dataset_name_new with 
     its exact length with the \0 termination */
  
  strncpy(dataset_name_new, dataset_name, *name_len);
  *(dataset_name_new + *name_len) = '\0';




  /* open the dataset */
  dataset = H5Dopen(*file_identifier, dataset_name_new); 

  dataspace = H5Dget_space(dataset);

  
  /* create the hyperslab -- this will differ on the different processors */
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_size);

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
                        &stride_1d, &count_1d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for generic_real_arr dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for generic_real_arr dataspace\n");
  }

  /* create the memory space */
  rank = 1;
  dimens_1d = *local_size;
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);

  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                H5P_DEFAULT, generic_arr);

  if (status < 0){
    printf("Error: Unable to read data for generic_real_arr\n");
    Driver_abortFlashC("Error: Unable to read data for generic_real_arr\n");
  }

  
  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  free(dataset_name_new);

}
