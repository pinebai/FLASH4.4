/* This file contains the routines that open and close the pnetcdf files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pnetcdf.h>
#include "mangle_names.h"
#include <mpi.h>
#include "constants.h"
/*#include "driverAPI.h"*/

/* define an info object to store MPI-IO information */
static MPI_Info FILE_INFO_TEMPLATE;


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_initialize_file)(int* file_identifier, 
                            char filename[MAX_STRING_LENGTH+1],
                            MPI_Fint* io_comm_,
                            int* io_outputSplitNum)
{
  MPI_Comm io_comm = MPI_Comm_f2c(*io_comm_);
  int ierr, status;

  /* make the filename pretty -- cutoff any trailing white space */
  int len = 0;
  char* string_index; 

  /* operate on a copy of the filename -- adding the null character for
     C messes with FORTRAN */

  char local_filename[MAX_STRING_LENGTH+1];


  string_index = filename;
  
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';
  

  /* ---------------------------------------------------------------------
      platform dependent code goes here -- the access template must be
      tuned for a particular filesystem blocksize.  some of these 
      numbers are guesses / experiments, others come from the file system
      documentation.

      The sieve_buf_size should be equal a multiple of the disk block size
     ---------------------------------------------------------------------- */

  /* create an MPI_INFO object -- on some platforms it is useful to
     pass some information onto the underlying MPI_File_open call */
  ierr = MPI_Info_create(&FILE_INFO_TEMPLATE);
  if(ierr < 0){
    printf("Error: Unable to create FILE_INFO_TEMPLATE\n");
    /*Driver_abort("Error: Unable to create FILE_INFO_TEMPLATE\n");*/
  }


#ifdef TFLOPS
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "access_style", "write_once");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "collective_buffering", "true");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", "1048576");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "4194304");
#endif

#ifdef SGI
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "access_style", "write_once");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "collective_buffering", "true");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", "1048576");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "4194304");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "sieve_buf_size", "4096");
#endif

#ifdef DEBUG_IO
  printf("set fapl to use MPI-IO, ierr = %d\n", (int) ierr);
#endif

  /* ----------------------------------------------------------------------
      end of platform dependent properties
     ---------------------------------------------------------------------- */

  /* create the file collectively */

  if(*io_outputSplitNum == 1){
    status = ncmpi_create(MPI_COMM_WORLD, local_filename, NC_CLOBBER | NC_64BIT_OFFSET, FILE_INFO_TEMPLATE, file_identifier);
  }else{
    status = ncmpi_create(io_comm, local_filename, NC_CLOBBER | NC_64BIT_OFFSET, FILE_INFO_TEMPLATE, file_identifier);
  }

  
  if(status < 0){
    printf("Error: ncmpi_create");
    Driver_abortFlashC("Error: ncmpi_create");
  }
  

#ifdef DEBUG_IO
        if (status != NC_NOERR)  handle_error(status);
#endif

#ifdef DEBUG_IO
  printf("openned the file, identifier = %d\n", (int) *file_identifier);
#endif

}



/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_open_file_for_read)(int *file_identifier, char filename[MAX_STRING_LENGTH+1])

{

  int len = 0;
  char* string_index;
  char local_filename[MAX_STRING_LENGTH+1];

  int status;

  MPI_Info_create(&FILE_INFO_TEMPLATE);
   

  string_index = filename;
 
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';
   
   status = ncmpi_open(MPI_COMM_WORLD, local_filename, NC_NOWRITE, FILE_INFO_TEMPLATE, file_identifier);
   if(status < 0){
     printf("Error: ncmpi_open");
     /*Driver_abort("Error: ncmpi_open");*/
   }
  

#ifdef DEBUG_IO
        if (status != NC_NOERR)  handle_error(status);
#endif


}     

void FTOC(io_ncmpi_close_file)(int* file_identifier)
{
  int status;

  /* close the file */
  status = ncmpi_close(*file_identifier);
  if(status < 0){
    printf("Error: ncmpi_close");
    /*Driver_abort("Error: ncmpi_close");*/
  }
  

#ifdef DEBUG_IO
        if (status != NC_NOERR)  handle_error(status);
#endif
  
  MPI_Info_free(&FILE_INFO_TEMPLATE);
  
}

