#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include <pnetcdf.h>
#include "mangle_names.h"
#include <mpi.h>
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);


/****if* source/IO/IOMain/pnetcdf/io_ncmpi_read_generic_iarr_
 *
 * NAME
 *
 *  io_ncmpi_read_generic_iarr_
 *
 *
 * SYNOPSIS
 *
 *  void io_ncmpi_read_generic_iarr_(file_identifier, generic_arr, local_size, total_size,
 *                        global_offset, dataset_name, dataset_name_len)
 *
 *  void io_ncmpi_read_generic_iarr_(hid_t*, int*, int*, int*, int*, char*, int*)
 *
 *
 * DESCRIPTION
 *
 *  Read in a single processors's portion of the data
 *  (generic_arr) to a FLASH parallel netcdf checkpoint file.  Given the ncmpi file handle 
 *  (file_identifier), the offset to start writing at (global_offset) read the
 *  local_size worth of data.  For this routine, often just one processor needs
 *  to read in the data.  This routine is typically called from 
 *  IO_readUserArray.  This routine reads in a 1 dimensional dataset.  
 *  Fortran multidimensional arrays will be read in as a 1d array.
 *
 *
 * ARGUMENTS
 *
 *  file_identifier       the ncmpi file handle for the current file
 *                        
 *
 *  generic_arr           user defined array to read to file
 *
 *  local_size            the local size of array to read
 *
 *  total_size            the total size(number of elements) of the array to read
 *                          
 *
 *  global_offset         the starting position of the current processor's
 *                        data (specified in blocks)
 * 
 *  dataset_name          name of the dataset to read out
 *
 *  dataset_name_len      length of the dataset name. makes things much easier
 *                        when trying to convert from c to fortran strings*  
 *  
 *
 ****/



/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_read_generic_iarr)(int* file_identifier,
				      int generic_arr[],
				      int* local_size,
				      int* total_size,
				      int* global_offset,
				      char* dataset_name, 
				      int* dataset_name_len)
{


  int ncid, status, varid;
  int dimid, rank;

  MPI_Offset start_1d, count_1d;

  char* dataset_name_new;
  char* dim_name_new;
  
  dataset_name_new = (char *) malloc((*dataset_name_len) + 1 * sizeof(char)); 
  
  /* copy the dataset_name into a c string dataset_name_new with 
     its exact length with the \0 termination */
  
  strncpy(dataset_name_new, dataset_name, *dataset_name_len);
  *(dataset_name_new + *dataset_name_len) = '\0';



  ncid = *file_identifier;
  

  /*defining a new dimension
    if the dimension is already defined the user could change this piece of code to say */
    
  status = ncmpi_inq_varid(ncid, dataset_name_new, &varid); 

    
  if (status < 0){
    printf("Error: Unable to get varid io_ncmpi_read_generic_arr\n");
    Driver_abortFlashC("Error: Unable to get varid io_ncmpi_read_generic_arr\n");
  }    



  start_1d = (MPI_Offset) (*global_offset);
  count_1d = (MPI_Offset) (*local_size);

  

  /* High Level Data Mode */
  status = ncmpi_get_vara_int_all(ncid, varid, &start_1d, &count_1d, generic_arr);



  if (status < 0){
    printf("Error: ncmpi_get_vara_int_all, generic_arr\n");
    Driver_abortFlashC("Error: ncmpi_get_vara_int_all, generic_arr\n");
    
  }

  free(dataset_name_new);
#ifdef DEBUG_IO
        if (status != NC_NOERR)  handle_error(status);
#endif
	
	


}
