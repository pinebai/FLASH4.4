/* This file contains the routines that open and close the HDF5 files */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
/*#include "driverAPI.h"*/

int HDF5_MODE = 0;		/* For serial HDF5 this should never become COLLECTIVE - KW */
/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/serial/io_h5init_file_
 * NAME
 *
 *  io_h5init_file_
 *
 *
 * SYNOPSIS
 *
 *  io_h5init_file_(file_identifier, filename)
 *
 *  io_h5init_file_(hid_t *, char [])
 *
 *
 * DESCRIPTION
 *
 *  Open an HDF5 file for output.  If successful in openning the file for 
 *  output, a positive file_identifier will be returned via the argument list.  
 *  This file_identifier is needed in all subsequent HDF5 calls that write 
 *  to the file.
 *
 *
 * NOTES
 *
 *  Requires HDF5 v. 1.4.0 or later.
 *
 ****/

void FTOC(io_h5init_file)(hid_t* file_identifier, char filename[MAX_STRING_LENGTH+1])
{


  /* make the filename pretty -- cutoff any trailing white space */
  int len = 0;
  char* string_index; 


  /* operate on a copy of the filename -- adding the null character for
     C messes with FORTRAN */

  char local_filename[MAX_STRING_LENGTH+1];

  string_index = filename;

  /* copy the filename into a buffer and make it null terminated  -- we 
     must be very careful here.  filename[] as it is passed in will not be
     null terminated (since Fortran doesn't believe in null), and the only
     way to find the end of the string is to search for the first occurance
     of whitespace (since we expect of filenames to be absent of any such 
     space).  Using a strcpy here can do very bad things, because there is 
     not \0 */


  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';


  /* open the file */
  *file_identifier = H5Fcreate(local_filename, H5F_ACC_TRUNC, 
                               H5P_DEFAULT, H5P_DEFAULT);

}


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/serial/io_h5open_file_for_read_
 *
 * NAME
 *
 *  io_h5open_file_for_read_
 *
 *
 * SYNOPSIS
 *
 *  io_h5open_file_for_read_(file_identifier, filename)
 *
 *  io_h5open_file_for_read_(hid_t *, char [])
 *
 *
 * DESCRIPTION
 *
 *  Open an existing HDF5 file for reading.  If successful, the 
 *  file_identifier will be positive when returned.
 *
 *
 * NOTES
 *
 *  This routine requires HDF5 v. 1.4.0 or later
 *
 ****/


void FTOC(io_h5open_file_for_read)(hid_t* file_identifier, char filename[MAX_STRING_LENGTH+1])
{


  /* make the filename pretty -- cutoff any trailing white space */
  int len = 0;
  char* string_index; 
  char local_filename[MAX_STRING_LENGTH+1];
    
  string_index = filename;

  while (*string_index != ' ') {
    local_filename[len] = filename[len];
    len++;
    string_index++;
  }

  *(local_filename+len) = '\0';

  /* open the file */
  *file_identifier = H5Fopen(local_filename, H5F_ACC_RDONLY, H5P_DEFAULT); 
    
  if (*file_identifier < 0) {
    printf("[H5_OPEN_FILE_FOR_READ] ERROR: opening file %s for input", local_filename);
    /*abort_flash("[H5_OPEN_FILE_FOR_READ] ERROR: opening file for input\n");*/
  }
}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/serial/io_h5close_file_
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
  /* close the file */
  H5Fclose(*file_identifier);
    
}










