/* This file contains the routines that open and close the HDF5 files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include <mpi.h>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>


/* define an info object to store MPI-IO information */
static MPI_Info FILE_INFO_TEMPLATE;
int HDF5_MODE = 0;
/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */


/****if* source/IO/IOMain/hdf5/parallel/h5_file_interface/io_h5set_xfer_mode_
 * NAME
 *
 *  io_h5set_xfer_mode_
 *
 *
 * SYNOPSIS
 *
 *  io_h5set_xfer_mode_(int *)
 *
 * DESCRIPTION
 *
 *  Sets HDF5_MODE to INDEPENDENT or COLLECTIVE.  This function is called
 *  from IO_Init.
 *
 ****/

void FTOC(io_h5set_xfer_mode)(const int * const pCollectiveHDF5)
{
  const int collectiveHDF5 = *pCollectiveHDF5;

  /*Set whether or not we are using collective mode or not.
    Done to allow us to set this with a runtime parameter*/
  HDF5_MODE = collectiveHDF5;
}


/****if* source/IO/IOMain/hdf5/parallel/h5_file_interface/io_h5init_file_
 * NAME
 *
 *  io_h5init_file_
 *
 *
 * SYNOPSIS
 * 
 *  io_h5init_file_(file_identifier, filename, io_comm, io_outputSplitNum, existing)
 *
 *  io_h5init_file_(hid_t *, char [], MPI_Comm* , int*, int*)
 *
 *
 * DESCRIPTION
 *
 *  Open an HDF5 file for output via parallel I/O (using the underly MPI-IO
 *  library).  Also set the buffer size and whether alignment to filesystem
 *  blocks is desired.  On some platforms, the MPI_Info object is used
 *  to pass hints to the underlying MPI-IO layer.
 *
 *  If successful in openning the file for output, a positive file_identifier
 *  will be returned via the argument list.  This file_identifier is needed
 *  in all subsequent HDF5 calls that write to the file.
 *
 *  io_comm and io_outputSplitNum are for IO file splitting.  This
 *  functionality is not yet fully implemented.  Intended to create
 *  multiple io communicators so users can break checkpoint files up
 *  into groups.
 *
 *  The existing argument is used to specify that this is an existing
 *  file that should be reopened so that it can be appended to.
 *  
 * NOTES
 * 
 *  Requires HDF5 v. 1.4.0 or later.
 *
 ****/

void FTOC(io_h5init_file)(hid_t* file_identifier, 
                    char filename[MAX_STRING_LENGTH+1],
                    MPI_Fint* io_comm_,
		    int* io_outputSplitNum,
		    int* existing)
{

  MPI_Comm io_comm = MPI_Comm_f2c(*io_comm_);
  int ierr;

  hid_t acc_template;

  /* make the filename pretty -- cutoff any trailing white space */
  int len = 0;
  char* string_index; 

  /* operate on a copy of the filename -- adding the null character for
     C messes with FORTRAN */
  char local_filename[MAX_STRING_LENGTH+1];

  string_index = filename;

  /* copy the filename into a buffer and make it null terminated */
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';


  /* set the file access template for parallel IO access */
  acc_template = H5Pcreate(H5P_FILE_ACCESS);

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

#ifdef IBM
  ierr = H5Pset_sieve_buf_size(acc_template, 262144); 
  ierr = H5Pset_alignment(acc_template, 524288, 262144);
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "IBM_largeblock_io", "true"); 
#endif

#ifdef TFLOPS
  ierr = H5Pset_sieve_buf_size(acc_template, 524288); 
  ierr = H5Pset_alignment(acc_template, 524288, 262144);

  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "access_style", "write_once");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "collective_buffering", "true");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", "1048576");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "4194304");
#endif

#ifdef CHIBA
  ierr = H5Pset_sieve_buf_size(acc_template, 524288); 
  ierr = H5Pset_alignment(acc_template, 524288, 262144);

  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "access_style", "write_once");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "collective_buffering", "true");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_block_size", "1048576");
  ierr = MPI_Info_set(FILE_INFO_TEMPLATE, "cb_buffer_size", "4194304");
#endif

  /* if we are splitting the output file ... */
  /*if(*io_outputSplitNum == 1){
    ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, FILE_INFO_TEMPLATE);
  }else{
  */ierr = H5Pset_fapl_mpio(acc_template, io_comm, FILE_INFO_TEMPLATE);
  /*}*/

#ifdef DEBUG_IO
  printf("set fapl to use MPI-IO, ierr = %d\n", (int) ierr);
#endif

  /* ----------------------------------------------------------------------
      end of platform dependent properties
     ---------------------------------------------------------------------- */


  /* ----------------------------------------------------------------------
      the jury is still out on dataset chunking, so group all of the 
      chunking settings under the CHUNK preprocessor variable
     ---------------------------------------------------------------------- */
#ifdef CHUNK
  ierr = H5Pset_cache(acc_template, 128, 512, (size_t) 1048576, 1.0);
#endif
  
  /* create the file collectively */
  if(*existing == 0) {
    *file_identifier = H5Fcreate(local_filename, H5F_ACC_TRUNC, 
				 H5P_DEFAULT, acc_template);
  } else {
    *file_identifier = H5Fopen(local_filename, H5F_ACC_RDWR, acc_template);
  }

#ifdef DEBUG_IO
  printf("openned the file, identifier = %d\n", (int) *file_identifier);
#endif

  /* release the file access template */
  ierr = H5Pclose(acc_template);

    
}



/****if* source/IO/IOMain/hdf5/parallel/h5_file_interface/io_h5open_file_for_read_
 *
 * NAME
 * 
 *  io_h5open_file_for_read_
 *
 *
 * SYNOPSIS
 *
 *  io_h5open_file_for_read_(file_identifier, filename, io_comm, io_outputSplitNum)
 *
 *  io_h5open_file_for_read_(hid_t *, char [], MPI_comm*, int*)
 *
 *
 * DESCRIPTION
 *
 *  Open an existing HDF5 file for reading using parallel I/O.  If 
 *  successful, the file_identifier will be positive when returned.
 *
 *  Last 2 args are not fully implemented.  For file splitting.
 *
 * NOTES
 *
 *  This routine requires HDF5 v. 1.4.0 or later
 *
 ****/ 

void FTOC(io_h5open_file_for_read)(hid_t* file_identifier, 
                           char filename[MAX_STRING_LENGTH+1],
                           MPI_Fint* io_comm_,
                           int* io_outputSplitNum)
{

  hid_t acc_template;
  MPI_Comm io_comm = MPI_Comm_f2c(*io_comm_);


  /* make the filename pretty -- cutoff any trailing white space */
  int len = 0;
  char* string_index; 
  char local_filename[MAX_STRING_LENGTH+1];

  int ierr;
  int FAIL = -1;

  /* create an MPI_INFO object -- on some platforms it is useful to
     pass some information onto the underlying MPI_File_open call */
  ierr = MPI_Info_create(&FILE_INFO_TEMPLATE);
    

  string_index = filename;
  
  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  acc_template = H5Pcreate(H5P_FILE_ACCESS);

 

  /* if we are splitting the output file ... */
  /*if(*io_outputSplitNum == 1){
     ierr = H5Pset_fapl_mpio(acc_template, MPI_COMM_WORLD, MPI_INFO_NULL);
  }else{
    */ierr = H5Pset_fapl_mpio(acc_template, io_comm, MPI_INFO_NULL);
  /*}*/


  MPI_Barrier(MPI_COMM_WORLD);

  *file_identifier = H5Fopen(local_filename,H5F_ACC_RDONLY,acc_template);

  ierr = H5Pclose(acc_template);
    
  if (*file_identifier < 0) {
    printf("Error opening file %s for input", local_filename);
    ierr = MPI_Abort(MPI_COMM_WORLD, FAIL);
  }

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/parallel/h5_file_interface/io_h5close_file_
 *
 * NAME
 *
 *  io_h5close_file_
 *
 * 
 * SYNOPSIS
 *
 *  io_h5close_file_(file_identifier)
 *
 *  io_h5close_file(hid_t *)
 * 
 *
 * DESCRIPTION
 *
 *  Close the HDF5 file pointed to by file_identifier.
 *
 ****/

void FTOC(io_h5close_file)(hid_t* file_identifier)
{
  int ierr;


  /* close the file */
  H5Fclose(*file_identifier);
  
  ierr = MPI_Info_free(&FILE_INFO_TEMPLATE);
  
}









