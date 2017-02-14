/****f* source/Driver/Driver_abortFlashC
*
* NAME
*
*  Driver_abortFlashC
*
* SYNOPSIS
*
*  Driver_abortFlashC(char* :: errorMessage)
*
* DESCRIPTION
*
*  C specific function that writes an error message 
*  to the logfile and aborts FLASH.
*  Attempts to shut down all processes (using MPI_Abort()).
*  
*  Use this function to abort the FLASH code from a C routine
*  
*
* ARGUMENTS
*
*  errorMessage:     A string to write to the logfile (presumably 
*                    indicating what went wrong).
*
* NOTES
*  for those familiar with FLASH2, we added this function to make it
*  easier to compile FLASH across different compilers.  It is possible
*  to call fortran90 routines from C functions, however an interface
*  is needed and some complications were found on various compilers.
*  Often users had to remove all the calls to abort_flash in the io 
*  C routines in order to get FLASH2 to run.  
*
*
***/


#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "constants.h"

int Driver_abortFlashC(char* message){


  int ierr;


  ierr = 0;
  
  return ierr;

}
  
  
