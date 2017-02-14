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


/* 
   This function writes a single unknown (passed from the checkpoint 
   routine).  Each variable has a corresponding varid in parallel netcdf.
   
   The dimensions of the unknowns array (nxb, nyb, nzb, local_blocks)
   are passed through as arguments.  We also store the maximum and minimum 
   values of the datasets.  
*/




/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_write_unknowns_sp)(int* file_identifier,
				      int* varid,
				      int* nxb,             /* num of zones to store in x */
				      int* nyb,             /* num of zones to store in y */
				      int* nzb,             /* num of zones to store in z */
				      float* unknowns,     /* [nvar][NZB][NYB][NXB] */
				      float* varMin, 
				      float* varMax,
				      int* local_blocks,
				      int* total_blocks,
				      int* global_offset)
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


  status = ncmpi_put_vara_float_all(ncid, *varid, start_4d, count_4d, unknowns);
  if (status < 0){
    printf("Error: ncmpi_put_vara_float_all, unknowns\n");
    Driver_abortFlashC("Error: ncmpi_put_vara_float_all, unknowns\n");
    
  }
  

  /* write out max and minimum values of the dataset */
  /*re enter define mode */
  status = ncmpi_redef(ncid);


  if (status < 0){
    printf("Error: write_unknowns_sp redef\n");
    Driver_abortFlashC("Error: write_unknowns_sp redef\n");
  }

  
  status = ncmpi_put_att_float(ncid, *varid, "minimum", NC_FLOAT, 1, varMin);
  if (status < 0){
    printf("Error: ncmpi_put_att_double minimum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, minimum\n");
  }
  
  status = ncmpi_put_att_float(ncid, *varid, "maximum", NC_FLOAT, 1, varMax);
  if (status < 0){
    printf("Error: ncmpi_put_att_double maximum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, maximum\n");
  }
  
  
  /* end define mode */
  status = ncmpi_enddef(ncid);

  if (status < 0){
    printf("Error: write_unknowns_sp enddef\n");
    Driver_abortFlashC("Error: write_unknowns_sp enddef\n");
  }

 
}

