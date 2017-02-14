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


/****if* source/IO/IOMain/pnetcdf/io_ncmpi_write_generic_iarr_
 *
 * NAME
 *
 *  io_ncmpi_write_generic_iarr_
 *
 *
 * SYNOPSIS
 *
 *  void io_ncmpi_write_generic_iarr_(file_identifier, generic_arr, local_size, total_size,
 *                        global_offset, dataset_name, dataset_name_len, dim_name, dim_name_len)
 *
 *  void io_ncmpi_write_generic_iarr_(hid_t*, int*, int*, int*, int*, char*, int*, char*, int*)
 *
 *
 * DESCRIPTION
 *
 *  Write out a single processors's portion of the data
 *  (generic_arr) to a FLASH parallel netcdf checkpoint file.  Given the ncmpi file handle 
 *  (file_identifier), the offset to start writing at (global_offset) write the
 *  local_size worth of data.  For this routine, often just one processor needs
 *  to write out the data.  This routine is typically called from 
 *  IO_writeUserArray.  This routine writes out a 1 dimensional dataset.  
 *  Fortran multidimensional arrays will be written as a 1d array.
 *
 *
 * ARGUMENTS
 *
 *  file_identifier       the ncmpi file handle for the current file
 *                        
 *
 *  generic_arr           user defined array to write to file
 *
 *  local_size            the local size of array to write
 *
 *  total_size            the total size(number of elements) of the array to write
 *                          
 *
 *  global_offset         the starting position of the current processor's
 *                        data (specified in blocks)
 * 
 *  dataset_name          name of the dataset to write out
 *
 *  dataset_name_len      length of the dataset name. makes things much easier
 *                        when trying to convert from c to fortran strings*  
 *  
 *  dim_name              name of the dimension to write out
 *
 *  dim_name_len          length of the dimension name. makes things much easier
 *                        when trying to convert from c to fortran strings
 *
 ****/


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_write_generic_iarr)(int* file_identifier,
				       int generic_arr[],
				       int* local_size,
				       int* total_size,
				       int* global_offset,
				       char* dataset_name, 
				       int* dataset_name_len,
				       char* dim_name, 
				       int* dim_name_len)
{


  int ncid, status, nbytes, dim_generic, nvars, varid;
  int dimid, rank;

  MPI_Offset start_1d, count_1d;

  char* dataset_name_new;
  char* dim_name_new;
  
  dataset_name_new = (char *) malloc((*dataset_name_len) + 1 * sizeof(char)); 
  
  /* copy the dataset_name into a c string dataset_name_new with 
     its exact length with the \0 termination */
  
  strncpy(dataset_name_new, dataset_name, *dataset_name_len);
  *(dataset_name_new + *dataset_name_len) = '\0';



  dim_name_new = (char *) malloc((*dim_name_len) + 1 * sizeof(char)); 
  
  /* copy the dataset_name into a c string dataset_name_new with 
     its exact length with the \0 termination */
  
  strncpy(dim_name_new, dim_name, *dim_name_len);
  *(dim_name_new + *dim_name_len) = '\0';



  ncid = *file_identifier;
  
  
  /*re enter define mode */
  status = ncmpi_redef(ncid);


  /*defining a new dimension
    if the dimension is already defined the user could change this piece of code to say
    
    status = ncmpi_inq_dimid(ncid, "dim_tot_blocks", &dim_tot_blocks); 
  */


  status = ncmpi_def_dim(ncid, dim_name_new, (MPI_Offset)(*total_size), &dim_generic);
    
  if (status < 0){
    printf("Error: Unable to define dim in write io_ncmpi_write_generic_arr\n");
    Driver_abortFlashC("Error: Unable to define dim in write io_ncmpi_write_generic_arr\n");
  }    


  status = ncmpi_inq_nvars(ncid, &nvars);
  varid = nvars;
    

  /* define var for generic array */
  rank = 1;
  
  status = ncmpi_def_var (ncid, dataset_name_new, NC_INT, rank, &dim_generic, &varid);
  
  if (status < 0){
    printf("Error: Unable to define var in io_ncmpi_write_generic_arr\n");
    Driver_abortFlashC("Error: Unable to define var io_ncmpi_write_generic_arr\n");
  }  

 
  
  status = ncmpi_enddef(ncid);
  
  if (status < 0){
    printf("Error io_ncmpi_write_generic_iarr: can not enddef\n");
    Driver_abortFlashC("Error io_ncmpi_write_generic_iarr: enddef\n");
  }
 






  start_1d = (MPI_Offset) (*global_offset);
  count_1d = (MPI_Offset) (*local_size);

  /* nbytes is size of int * the number of dimensions, in the case of 
     generic_arr there is only one dimension, num blocks
  */
  nbytes = (int)sizeof(int);
  

  /* High Level Data Mode */
  status = ncmpi_put_vara_int_all(ncid, varid, &start_1d, &count_1d, generic_arr);

  /* flexible Data Mode Interface way 
  status = ncmpi_put_vara_all(ncid, *varid, &start_1d, &count_1d, 
                        (const void*)generic_arr, nbytes, MPI_INT);*/


  if (status < 0){
    printf("Error: ncmpi_put_vara_int_all, generic_arr\n");
    Driver_abortFlashC("Error: ncmpi_put_vara_int_all, generic_arr\n");
    
  }

  free(dataset_name_new);

#ifdef DEBUG_IO
        if (status != NC_NOERR)  handle_error(status);
#endif

 /* end define mode */


}
