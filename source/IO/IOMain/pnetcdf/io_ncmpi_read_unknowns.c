#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <pnetcdf.h>
#include "mangle_names.h"
#include <mpi.h>
#include <assert.h>
#include "Flash.h"
#include "constants.h"



/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/* 
   This function reads in a single unknown (passed from the checkpoint 
   routine).  Each variable has a corresponding varid in parallel netcdf.

   
   The dimensions of the unknowns array (nxb, nyb, nzb, local_blocks)
   are passed through as arguments.  
*/

void FTOC(io_ncmpi_read_unknowns)(int* file_identifier,
                         int* varid,
                         int* nxb,             /* num of zones to store in x */
                         int* nyb,             /* num of zones to store in y */
                         int* nzb,             /* num of zones to store in z */
                         double* unknowns,     /* [mblks][NZB][NYB][NXB] */
                         int* local_blocks,  /* local num blocks to read */
                         int* total_blocks,   /* not used ... */
                         int* global_offset)   /*offset for current proc to start reading */
{
  int ncid, status;
  
  MPI_Offset start_4d[4], count_4d[4];
  
  ncid = *file_identifier;
  
  /* create the hyperslab -- this will differ on the different processors */
  start_4d[0] = (MPI_Offset) (*global_offset);
  start_4d[1] = 0;
  start_4d[2] = 0;
  start_4d[3] = 0;
  
  count_4d[0] = (MPI_Offset) (*local_blocks);
  count_4d[1] = *nzb;
  count_4d[2] = *nyb;
  count_4d[3] = *nxb;
  
  status = ncmpi_get_vara_double_all(ncid, *varid, start_4d, count_4d, unknowns);
  if (status < 0){
    printf("Error: ncmpi_get_vara_double_all, unknowns\n");
    Driver_abortFlashC("Error: ncmpi_get_vara_double_all, unknowns\n");
  }

#ifdef DEBUG_IO
  if (status != NC_NOERR)  handle_error(status);
#endif

}

