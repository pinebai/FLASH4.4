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

void handle_error(int status) {
  /*fprintf(stderr, "%s\n", ncmpi_strerror(status));*/
}





/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */


/* 
   This function writes a single unknown (passed from the checkpoint 
   routine).  Each variable has a corresponding varid in parallel netcdf.
   
   The dimensions of the unknowns array (nxb, nyb, nzb, local_blocks)
   are passed through as arguments.  We also store the maximum and minimum 
   values of the datasets.  
*/


void FTOC(io_ncmpi_write_unknowns)(int* file_identifier,
				   int* varid,
				   int* nxb,             /* num of zones to store in x */
				   int* nyb,             /* num of zones to store in y */
				   int* nzb,             /* num of zones to store in z */
				   double* unknowns,     /* [nvar][NZB][NYB][NXB] */
				   double* varMin,
				   double* varMax,
				   int* local_blocks,
				   int* total_blocks,
				   int* global_offset)
{
  int ncid, status;
  
  MPI_Offset start_4d[4], count_4d[4];
  
  ncid = *file_identifier;
  
  /*
    printf("local_blocks = %d\n", *local_blocks);
    printf("total_blocks = %d\n", *total_blocks);
    printf("global_offset = %d\n", *global_offset);
    printf("varid = %d\n", *varid);
  */

  /* create the hyperslab -- this will differ on the different processors */
  start_4d[0] = (MPI_Offset) (*global_offset);
  start_4d[1] = 0;
  start_4d[2] = 0;
  start_4d[3] = 0;
  
  count_4d[0] = (MPI_Offset) (*local_blocks);
  count_4d[1] = *nzb;
  count_4d[2] = *nyb;
  count_4d[3] = *nxb;


  status = ncmpi_put_vara_double_all(ncid, *varid, start_4d, count_4d, unknowns);
  if (status < 0){
    printf("Error: ncmpi_put_vara_double_all, unknowns\n");
    Driver_abortFlashC("Error: ncmpi_put_vara_double_all, unknowns\n");
    
  }


  /* write out max and minimum values of the dataset */
  /*re enter define mode */
  status = ncmpi_redef(ncid);


  if (status < 0){
    printf("Error: write_unknowns redef\n");
    Driver_abortFlashC("Error: write_unknowns redef\n");
  }


  status = ncmpi_put_att_double(ncid, *varid, "minimum", NC_DOUBLE, 1, varMin);
  if (status < 0){
    printf("Error: ncmpi_put_att_double minimum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, minimum\n");
  }

  status = ncmpi_put_att_double(ncid, *varid, "maximum", NC_DOUBLE, 1, varMax);
  if (status < 0){
    printf("Error: ncmpi_put_att_double maximum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, maximum\n");
  }


  /* end define mode */
  status = ncmpi_enddef(ncid);

  if (status < 0){
    printf("Error: write_unknowns enddef\n");
    Driver_abortFlashC("Error: write_unknowns enddef\n");
  }
  
}

