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

void FTOC(io_ncmpi_write_part_dims)(int* file_identifier,
                            int* varid,
                            int* total_blocks,
                            int* npart_props, /* number of particle properties */
                            int* totalparticles,
                            int* particlesToCheckpoint,
                            char propNames[][OUTPUT_PROP_LENGTH])

{
  int ncid, status, dim_particles, dim_npart_props, dim_mdim, rank, dim_tot_blocks;
  int dimids[2]; 
  MPI_Offset start_2d[2], count_2d[2], stride_2d[2];
  int varidpart, string_size, i;
  char tmpstr[OUTPUT_PROP_LENGTH+1];
  char mkeys[40];
  char str[4];
  char *p;



  ncid = *file_identifier;

  /*if we are writing particles to checkpoint then*/
  /*re enter define mode */
  if(*particlesToCheckpoint == 1){
    status = ncmpi_redef(ncid);

    if(*totalparticles > 0){
      /*get the dimid because it has already been defined in io_ncmpi_write_header */
      status = ncmpi_inq_dimid(ncid, "dim_tot_blocks", &dim_tot_blocks);
    }
      
  }else{

  
    if(*totalparticles > 0) {
      status = ncmpi_def_dim(ncid, "dim_tot_blocks", (MPI_Offset)(*total_blocks), &dim_tot_blocks);
      if (status < 0){
      printf("Error: Unable to define dim_tot_blocks in write part dims\n");
      Driver_abortFlashC("Error: Unable to define dim_tot_blocks in write part dims\n");
      }    
    }
  }
   
  
  /* define the dimensions */ 
  if(*totalparticles > 0) {
    status = ncmpi_def_dim(ncid, "dim_particles", (MPI_Offset)(*totalparticles), &dim_particles);
    if (status < 0){
      printf("Error: Unable to define dim_particles\n");
      Driver_abortFlashC("Error: Unable to define dim_particles\n");
    }    
    
    status = ncmpi_def_dim(ncid, "dim_npart_props", (MPI_Offset)(*npart_props), &dim_npart_props);
    if (status < 0){
      printf("Error: Unable to define dim_npart_props\n");
      Driver_abortFlashC("Error: Unable to define dim_npart_props\n");
    }    


    /* define var for localnp */
    rank = 1;
    dimids[0] = dim_tot_blocks;
    
    status = ncmpi_def_var (ncid, "localnp", NC_INT, rank, dimids, varid);
    
    if (status < 0){
      printf("Error: Unable to define local np var\n");
      Driver_abortFlashC("Error: Unable to define localnp var\n");
    }  



    /* define var for particles */
    rank = 2;
    dimids[0] = dim_particles;
    dimids[1] = dim_npart_props;
    
    varidpart = *varid + 1;
    
    status = ncmpi_def_var (ncid, "particles", NC_DOUBLE, rank, dimids, &varidpart);
    
    if (status < 0){
      printf("Error: Unable to define particles\n");
      Driver_abortFlashC("Error: Unable to define particles\n");
    }  
    
    
    
    
    /*write out the particle attribute labels */
    
    for(i=0; i<NPART_PROPS; i++) {
      sprintf(str, "%d", i);
      strcpy(mkeys, "particle_props_name_");
      strncpy(tmpstr,propNames[i], OUTPUT_PROP_LENGTH);
      strcat(mkeys,str);
      p = strtok(tmpstr, "    ");
      if (p == NULL) {
      p = tmpstr;
      p[0] = '\0';
      }
      string_size=strlen(p);
      
      status = ncmpi_put_att_text(ncid, NC_GLOBAL, mkeys, string_size, p);
      if (status < 0){
      printf("Error: Unable to write particle_props_name\n");
      Driver_abortFlashC("Error: Unable to write particle_props_name\n");
      }
    }
    
  }
      




  
  /* end define mode */

    status = ncmpi_enddef(ncid);
  
    if (status < 0){
      printf("Error io_ncmpi_write_part_dims: can not enddef\n");
      Driver_abortFlashC("Error io_ncmpi_write_part_dims: enddef\n");
    }

}
