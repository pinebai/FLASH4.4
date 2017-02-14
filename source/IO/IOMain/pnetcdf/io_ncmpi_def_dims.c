
/* if the following flag is defined, status checking on the writes will
   be performed, and the results will be output to stdout */

/* #define DEBUG_IO */

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

/*define dimensions all at once if possible.  pnetcdf has two modes "define mode" and 
  "write mode."  You could go back and forth between the 2 modes for each dataset to 
   write, but since we know most of the dimensions by the time this routine is called,
   we just write them all now for efficiency.  

   In contrast, for a few datasets we define the dimensions at the same time as
   writing the data.  We do this for particles for example, because they are not
   always included in a simulation
*/


void FTOC(io_ncmpi_def_dims)(int MyPE,                         
                   int nvar_out, 
                   int nzones_block[3],   
                   char unk_labels[][4],
                   int ncid,
                   int total_blocks,
                   int ngid)
     
{

  int i,k;
  int status = 0;
  int rank;
  int dim_tot_blocks, dim_nxb, dim_nyb, dim_nzb, dim_NDIM, dim_NGID, dim_2, dim_1;
  int dim_MAX_STR_LEN, dim_particles, dim_MDIM;
  int dimids[4]; 
  int two = 2;
  int one = 1;
  char *p;
  /* DEV: keep this until we have particles implemented */
  int nParticles = -1;


    
  status = ncmpi_def_dim(ncid, "dim_tot_blocks", (MPI_Offset)(total_blocks), &dim_tot_blocks);
  if (status < 0){
    printf("Error: Unable to define dim_tot_blocks\n");
    Driver_abortFlashC("Error: Unable to define dim_tot_blocks\n");
  }
 
  status = ncmpi_def_dim(ncid, "dim_nxb", (MPI_Offset)nzones_block[0], &dim_nxb);
  if (status < 0){
    printf("Error: Unable to define dim_nxb\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nyb", (MPI_Offset)nzones_block[1], &dim_nyb);
  if (status < 0){
    printf("Error: Unable to define dim_nyb\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nzb", (MPI_Offset)nzones_block[2], &dim_nzb);
  if (status < 0){
    printf("Error: Unable to define dim_nzb\n");
  }

  status = ncmpi_def_dim(ncid, "dim_NGID", (MPI_Offset)ngid, &dim_NGID);
  if (status < 0){
    printf("Error: Unable to define dim_NGID\n");
    Driver_abortFlashC("Error: Unable to define dim_NGID\n");
  }

  status = ncmpi_def_dim(ncid, "dim_NDIM", (MPI_Offset)NDIM, &dim_NDIM);
  if (status < 0){
    printf("Error: Unable to define dim_NDIM\n");
    Driver_abortFlashC("Error: Unable to define dim_NDIM\n");
  }

  status = ncmpi_def_dim(ncid, "dim_2", (MPI_Offset)two, &dim_2);
  if (status < 0){
    printf("Error: Unable to define dim_2\n");
    Driver_abortFlashC("Error: Unable to define dim_2\n");
  }


  status = ncmpi_def_dim(ncid, "dim_1", (MPI_Offset)one, &dim_1);
  if (status < 0){
    printf("Error: Unable to define dim_1\n");
    Driver_abortFlashC("Error: Unable to define dim_1\n");
  }


  status = ncmpi_def_dim(ncid, "dim_MAX_STR_LEN", (MPI_Offset)MAX_STRING_LENGTH, &dim_MAX_STR_LEN);
  if (status < 0){
    printf("Error: Unable to define dim_MAX_STR_LEN\n");
    Driver_abortFlashC("Error: Unable to define dim_MAX_STR_LEN\n");
  }

  status = ncmpi_def_dim(ncid, "dim_MDIM", (MPI_Offset)MDIM, &dim_MDIM);
  if(status < 0){
    printf("Error: Unable to define dim_MDIM\n");
    Driver_abortFlashC("Error: Unable to define dim_MDIM\n");
  }

}

























