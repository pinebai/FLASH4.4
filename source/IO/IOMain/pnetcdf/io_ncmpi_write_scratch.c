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





/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_ncmpi_write_scratch)(int* file_identifier,
				  int* nxb,             /* num of zones to store in x */
				  int* nyb,             /* num of zones to store in y */
				  int* nzb,             /* num of zones to store in z */
				  double* scratch,     /* [nvar][NZB][NYB][NXB] */
				  double* varMin,
				  double* varMax,
				  char record_label[5], /* add char-null termination */
				  int* local_blocks,
				  int* total_blocks,
				  int* global_offset)
{

  int ncid, status, rank, nvars, varid;
  int dim_tot_blocks, dim_nxb, dim_nyb, dim_nzb;
  int dimids[4]; 
  MPI_Offset start_4d[4], count_4d[4];


  char record_label_new[5];


  ncid = *file_identifier;

  /* 
     the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination 
  */
  strncpy(record_label_new, record_label,4);
  *(record_label_new + 4) = '\0';
  
  
  /*re enter define mode */
  status = ncmpi_redef(ncid);


  status = ncmpi_inq_nvars(ncid, &nvars);
  varid = nvars;
    

  /* define var for generic array */
  rank = 4;

  /*getting the already defined dimension */
    
  status = ncmpi_inq_dimid(ncid, "dim_tot_blocks", &dim_tot_blocks);
  if(status < 0){
    Driver_abortFlashC("Error: unable to get dim in io_ncmpi_write_scratch");
  }
  status = ncmpi_inq_dimid(ncid, "dim_nxb", &dim_nxb); 
  if(status < 0){
    Driver_abortFlashC("Error: unable to get dim in io_ncmpi_write_scratch");
  }  
  status = ncmpi_inq_dimid(ncid, "dim_nyb", &dim_nyb); 
  if(status < 0){
    Driver_abortFlashC("Error: unable to get dim in io_ncmpi_write_scratch");
  }  
  status = ncmpi_inq_dimid(ncid, "dim_nzb", &dim_nzb); 
  if(status < 0){
    Driver_abortFlashC("Error: unable to get dim in io_ncmpi_write_scratch");
  }

  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_nzb;
  dimids[2] = dim_nyb;
  dimids[3] = dim_nxb;
  
  status = ncmpi_def_var (ncid, record_label_new, NC_DOUBLE, rank, dimids, &varid);
  
  if (status < 0){
    printf("Error: Unable to define var in io_ncmpi_write_scratch\n");
    Driver_abortFlashC("Error: Unable to define var io_ncmpi_write_scratch\n");
  }  

  status = ncmpi_put_att_double(ncid, varid, "minimum", NC_DOUBLE, 1, varMin);
  if (status < 0){
    printf("Error: ncmpi_put_att_double minimum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, minimum\n");
  }
  
  status = ncmpi_put_att_double(ncid, varid, "maximum", NC_DOUBLE, 1, varMax);
  if (status < 0){
    printf("Error: ncmpi_put_att_double maximum");
    Driver_abortFlashC("Error: ncmpi_put_att_double, maximum\n");
  }

  
  status = ncmpi_enddef(ncid);
  
  if (status < 0){
    printf("Error io_ncmpi_write_scratch: can not enddef\n");
    Driver_abortFlashC("Error io_ncmpi_write_scratch: enddef\n");
  }

  /* create the hyperslab -- this will differ on the different processors */
  start_4d[0] = (MPI_Offset) (*global_offset);
  start_4d[1] = 0;
  start_4d[2] = 0;
  start_4d[3] = 0;
  
  count_4d[0] = (MPI_Offset) (*local_blocks);
  count_4d[1] = *nzb;
  count_4d[2] = *nyb;
  count_4d[3] = *nxb;


  status = ncmpi_put_vara_double_all(ncid, varid, start_4d, count_4d, scratch);
  if (status < 0){
    printf("Error: ncmpi_put_vara_double_all, scratch\n");
    Driver_abortFlashC("Error: ncmpi_put_vara_double_all, scratch\n");
    
  }
  
#ifdef DEBUG_IO
  if (status != NC_NOERR)  handle_error(status);
#endif
  
}

