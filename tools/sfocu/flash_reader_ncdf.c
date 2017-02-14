#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef NO_NCDF

#include "mpi.h"
#include "pnetcdf.h"
#include "flash_reader.h"
#include "options.h"

void FR_GetNumIntPartProps_NCDF(int *file_identifier,
				 int *numIntPartProps) {
  /* - maybe we'll need these variables again someday.
       Right now there are no integer particle properties
  MPI_Offset dimlen;
  int status;
  int dimid;
  */

  *numIntPartProps = 0;

}

void FR_GetNumRealPartProps_NCDF(int *file_identifier,
			     int *numRealPartProps) {
  MPI_Offset dimlen;
  int status;
  int dimid;

  status = ncmpi_inq_dimid(*file_identifier, "dim_npart_props", &dimid);

  if (status < 0){
    printf("Error: ncmpi_inq_dimid: dim_npart_props\n");
  }
  
  status = ncmpi_inq_dimlen(*file_identifier, dimid, &dimlen);
  if (status < 0){
    printf("Error: ncmpi_inq_dimlen: dim_npart_props\n");
  }
  *numRealPartProps = dimlen;
}


void handle_error(int status) {
  fprintf(stderr, "%s\n", ncmpi_strerror(status));
}

FR_ParticlesAllProps *FR_GetParticlesAllProps_NCDF(FR_File *file, int startingAt, 
						   int atMost, int *numberGotten) {
  MPI_Offset att_len;
  MPI_Offset start_2d[2], count_2d[2];
  int num_real_props;
  int varid;
  int status;
  int to_read;
  FR_ParticlesAllProps *particles;
  char att_name[50];
  int i;
  
  if (file->totalparticles > 0) {
    if (atMost < (file->totalparticles - startingAt)) {
      to_read = atMost;
    } else {
      to_read = file->totalparticles - startingAt;
    }

    particles = (FR_ParticlesAllProps *) malloc(sizeof(FR_ParticlesAllProps));

    particles->numIntProps = 0;

    FR_GetNumRealPartProps_NCDF((int*)file->handle, &num_real_props);
    particles->numRealProps = num_real_props;

    particles->realProps = (double *) malloc(sizeof(double) * to_read * particles->numRealProps);
    particles->intProps = (int *) malloc(sizeof(int) * to_read * particles->numIntProps);
    particles->realPropsNames = (char **) malloc(sizeof(char *) * particles->numRealProps);
    particles->intPropsNames  = (char **) malloc(sizeof(char *) * particles->numIntProps);

    for (i = 0; i < particles->numRealProps; i++) {
      particles->realPropsNames[i] = (char *) malloc(sizeof(char) * (FR_PART_PROP_STRING_SIZE+1));
      sprintf(att_name, "particle_props_name_%d", i);
      ncmpi_inq_attlen(*(int*)file->handle, NC_GLOBAL, att_name, &att_len);
      status = ncmpi_get_att_text(*(int*)file->handle, NC_GLOBAL, att_name, particles->realPropsNames[i]);
      particles->realPropsNames[i][att_len] = '\0';
      if (status < 0){
	printf("Error: ncmpi_get_att_int, p\n");
      }
    }

    start_2d[0] = (MPI_Offset)(startingAt);
    start_2d[1] = 0;
    
    count_2d[0] = (MPI_Offset)(to_read);
    count_2d[1] = (MPI_Offset)(particles->numRealProps);

    status = ncmpi_inq_varid(*(int*)file->handle, "particles", &varid);
    if (status < 0){
      printf("Error: ncmpi_inq_varid, real_props_id\n");
    }

    status = ncmpi_get_vara_double_all(*(int*)file->handle, varid, start_2d, count_2d, particles->realProps);
    if (status < 0){
      printf("Error: ncmpi_get_vara_real_all, particles->realProps\n");
    }

    *numberGotten = to_read;
    particles->numParticles = to_read;
  } else {
    *numberGotten = 0;
    particles = NULL;
  }
  return particles;
}


void FR_GetNumParticles_NCDF(int *file_identifier,
			     int *numParticles)
{

  /*
   * FR_FetNumParticles_NCDF: look at the particles entry and get the size
   *                of the dataset based on the dimensions of the array
   *
   * input:  file handle          (integer)
   *
   * output: number of particles  (integer)
   *
   * return value: -1 indicates an error getting the dimension
   *
   */

  MPI_Offset dimlen;
  int status;
  int dimid;

  *numParticles = 0;

  status = ncmpi_inq_dimid(*file_identifier, "dim_particles", &dimid);
  if (!(status < 0)){
  
    status = ncmpi_inq_dimlen(*file_identifier, dimid, &dimlen);
    if (status < 0){
      printf("Error: ncmpi_inq_dimlen: dim_particles\n");
    }

    *numParticles = dimlen;
  }
}




FR_File *FR_open_NCDF(char *filename, options_t opts){
  FR_File *out;
  int varid;
  int i, j;
  int status;
  char temp_varnames[FR_MAXVARS*FR_VAR_STRING_SIZE];
  int handle;
  int num_unk_names;
  char varname[FR_VAR_STRING_SIZE];
  char var_att_name[20];
  char att_name[50];
  MPI_Offset att_len;
  MPI_Offset start_1d, count_1d;
  MPI_Offset start_2d[2], stride_2d[2], count_2d[2];
  MPI_Offset start_3d[3], stride_3d[3], count_3d[3];

  MPI_Offset metadata_dim_size;
  int dimid;



  /* Open file */
  out = (FR_File *) malloc(sizeof(FR_File));
  out->handle = malloc(sizeof(int));
  out->format = FR_NCDF;
  strcpy(out->filename, filename);

  /*FIXME: we are not using the vartype in netcdf yet, this is a stopgap.*/
  for(i=0;i<FR_MAXVARS;++i){
   out->vartypes[i] = UNK;
  }

  status = ncmpi_open(MPI_COMM_WORLD, filename, NC_CLOBBER, MPI_INFO_NULL, &handle); 
  if(status < 0){
    printf("Error: ncmpi_open");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }
  *((int *)out->handle) = handle;

  /* total blocks */
  /* The 1 stands for an attribute of the 1st (counting from 0) 
  variable in this case the scalars variable, look in pnetcdf_write.c for
  more info */
  status = ncmpi_get_att_int(handle, 1, "globalnumblocks", &(out->nblocks)); 
  if (status < 0){
    printf("Error: ncmpi_get_att_int, globalNumBlocks\n");
    return NULL;
  }

  /* nxb, nyb, nzb */
  status = ncmpi_get_att_int(handle, 1, "nxb", &(out->ncells_vec[0])); 
  if (status < 0){
    printf("Error: ncmpi_get_att_int, nxb\n");
    return NULL;
  }

  status = ncmpi_get_att_int(handle, 1, "nyb", &(out->ncells_vec[1])); 
  if (status < 0){
    printf("Error: ncmpi_get_att_int, nyb\n");
    return NULL;
  }

  status = ncmpi_get_att_int(handle, 1, "nzb", &(out->ncells_vec[2])); 
  if (status < 0){
    printf("Error: ncmpi_get_att_int, nzb\n");
    return NULL;
  }

  /* use ncells_vec to get dim */
  if (out->ncells_vec[2] != 1)
    out->dim = 3;
  else if (out->ncells_vec[1] != 1)
    out->dim = 2;
  else
    out->dim = 1;

  /*now that we have the dimensionality, the size is either
    the current out->dim, or if dim_MDIM present, MDIM (3) */
  status = ncmpi_inq_dimid(handle, "dim_MDIM", &dimid);
  if(status != NC_NOERR){
    ncmpi_inq_dimid(handle, "dim_NDIM", &dimid);
  }
  
  ncmpi_inq_dimlen(handle, dimid, &metadata_dim_size);
  
  out->ncells = out->ncells_vec[0] *
                out->ncells_vec[1] *
                out->ncells_vec[2] ;
  

  out->nodetype = (int *) malloc(sizeof(int)*out->nblocks);  
  status = ncmpi_inq_varid(handle, "nodetype", &varid);
  if (status < 0){
    printf("Error: ncmpi_inq_varid, nodetype\n");
    return NULL;
  }
  
  start_1d = (MPI_Offset) 0;
  count_1d = (MPI_Offset) out->nblocks;

  status = ncmpi_get_vara_int_all(handle, varid, &start_1d, &count_1d, out->nodetype);
  if (status < 0){
    printf("Error: ncmpi_get_vara_int_all, nodetype\n");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }
  
  /*  for (i = 0; i < out->nblocks; i++) {
    printf("nodetype[%d]: %d\n", i, out->nodetype[i]);
  } 
  */
  out->lref = (int *) malloc(sizeof(int)*out->nblocks);  
  status = ncmpi_inq_varid(handle, "lrefine", &varid);
  if (status < 0){
    printf("Error: ncmpi_inq_varid, lrefine\n");
    return NULL;
  }
  
  start_1d = (MPI_Offset) 0;
  count_1d = (MPI_Offset) out->nblocks;

  status = ncmpi_get_vara_int_all(handle, varid, &start_1d, &count_1d, out->lref);
  if (status < 0){
    printf("Error: ncmpi_get_vara_int_all, lrefine\n");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }
  /*
  for (i = 0; i < out->nblocks; i++) {
    printf("lrefine[%d]: %d\n", i, out->lref[i]);
    } */



  out->coord = (double *) malloc(sizeof(double)*out->dim*out->nblocks);  
  status = ncmpi_inq_varid(handle, "coordinates", &varid);
  if (status < 0){
    printf("Error: ncmpi_inq_varid, coordinates\n");
    return NULL;
  }
  
  
  start_2d[0] = (MPI_Offset) 0;
  start_2d[1] = (MPI_Offset) 0;

  stride_2d[0] = (MPI_Offset) 1;
  stride_2d[1] = (MPI_Offset) 1;

  count_2d[0] = (MPI_Offset) out->nblocks;
  count_2d[1] = (MPI_Offset) out->dim;

  status = ncmpi_get_vars_double_all(handle, varid, start_2d, count_2d, stride_2d, out->coord);
  if (status < 0){
    printf("Error: ncmpi_get_vars_double_all, coord\n");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }

  /* for (i = 0; i < out->dim*out->nblocks; i++) {
    printf("coord[%d]: %f\n", i, out->coord[i]);
    }
  */

  out->size = (double *) malloc(sizeof(double)*out->dim*out->nblocks);  
  status = ncmpi_inq_varid(handle, "blocksize", &varid);
  if (status < 0){
    printf("Error: ncmpi_inq_varid, blocksize\n");
    return NULL;
  }
  
  start_2d[0] = (MPI_Offset) 0;
  start_2d[1] = (MPI_Offset) 0;

  stride_2d[0] = (MPI_Offset) 1;
  stride_2d[1] = (MPI_Offset) 1;

  count_2d[0] = (MPI_Offset) out->nblocks;
  count_2d[1] = (MPI_Offset) out->dim;

  status = ncmpi_get_vars_double_all(handle, varid, start_2d, count_2d, stride_2d, out->size);
  if (status < 0){
    printf("Error: ncmpi_get_vars_double_all, blocksize\n");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }
  /*
  for (i = 0; i < out->dim*out->nblocks; i++) {
    printf("size[%d]: %f\n", i, out->size[i]);
  } 
  */
  out->bbox = (double *) malloc(sizeof(double)*out->dim*out->nblocks*2);  
  status = ncmpi_inq_varid(handle, "bndbox", &varid);
  if (status < 0){
    printf("Error: ncmpi_inq_varid, bndbox\n");
    return NULL;
  }
  
  start_3d[0] = (MPI_Offset) 0;
  start_3d[1] = (MPI_Offset) 0;
  start_3d[2] = (MPI_Offset) 0;
 
  stride_3d[0] = (MPI_Offset) 1;
  stride_3d[1] = (MPI_Offset) 1;
  stride_3d[2] = (MPI_Offset) 1;

  count_3d[0] = (MPI_Offset) out->nblocks;
  count_3d[1] = (MPI_Offset) out->dim;
  count_3d[2] = (MPI_Offset) 2;

  status = ncmpi_get_vars_double_all(handle, varid, start_3d, count_3d,stride_3d,out->bbox);
  if (status < 0){
    printf("Error: ncmpi_get_vars_double_all, bndbox\n");
    if (status != NC_NOERR)  handle_error(status);
    return NULL;
  }

  /* number of variables and variable names */
  for(i=0; i<(FR_MAXVARS*FR_VAR_STRING_SIZE); i++)
    temp_varnames[i]='\0';

  /* unknown names */
  status = ncmpi_get_att_int(handle, NC_GLOBAL, "unknown_names", &num_unk_names); 
  if (status < 0){
    printf("Error: ncmpi_get_att_int, num_unk_names\n");
    return NULL;
  }

  out->nvar = num_unk_names;
  for (i = 0; i < num_unk_names; i++) {
    sprintf(var_att_name, "unknown_%d", i);
    status = ncmpi_get_att_text(handle, NC_GLOBAL, var_att_name, varname); 
    if (status < 0){
      printf("Error: ncmpi_get_att_text, num_unk_names\n");
      return NULL;
    }
    /* ncdf can't have ' 's in its var names, so 
       we use '_'s.  change them back to ' 's */
    for (j = 0; j < FR_VAR_STRING_SIZE; j++) {
      if (varname[j] == '_') {
	varname[j] = ' ';
      }
    }
    strncpy(out->varnames[i], varname, FR_VAR_STRING_SIZE);
    out->varnames[i][FR_VAR_STRING_SIZE] = '\0';
  }

  out->totalparticles = 0;
  out->numRealPartProps = 0;
  out->numIntPartProps = 0;
  FR_GetNumParticles_NCDF((int*)out->handle, &(out->totalparticles));
  if (out->totalparticles > 0) {
    FR_GetNumRealPartProps_NCDF((int*)out->handle, &(out->numRealPartProps));
    FR_GetNumIntPartProps_NCDF((int*)out->handle, &(out->numIntPartProps));

    for (i = 0; i < out->numIntPartProps; i++) {
      sprintf(att_name, "int_particle_props_name_%d", i);
      ncmpi_inq_attlen(*(int*)out->handle, NC_GLOBAL, att_name, &att_len);
      status = ncmpi_get_att_text(*(int*)out->handle, NC_GLOBAL, att_name, out->intPartPropNames[i]);
      out->intPartPropNames[i][att_len] = '\0';
      if (status < 0){
	printf("Error getting int part prop names\n");
      }
    }
    
    for (i = 0; i < out->numRealPartProps; i++) {
      sprintf(att_name, "particle_props_name_%d", i);
      ncmpi_inq_attlen(*(int*)out->handle, NC_GLOBAL, att_name, &att_len);
      status = ncmpi_get_att_text(*(int*)out->handle, NC_GLOBAL, att_name, out->realPartPropNames[i]);
      out->realPartPropNames[i][att_len] = '\0';
      if (status < 0){
	printf("Error getting real part prop names\n");
      }
    }
  }
  return out;  
}


FR_Block *FR_GetBlock_NCDF(FR_File *file, char *var, int block_no){
  FR_Block *out;
  int status;
  int varid;
  char mangled_name[FR_VAR_STRING_SIZE+1];
  int i;

  MPI_Offset start_4d[4], count_4d[4];

  strcpy(mangled_name, var);

  for (i = 0; i < FR_VAR_STRING_SIZE; i++) {
    if (mangled_name[i] == ' ') {
      mangled_name[i] = '_';
    }
  }

  out = (FR_Block *) malloc(sizeof(FR_Block));
  out->data = (double *) malloc(sizeof(double)*file->ncells);

	for(i=0; i < FR_MDIM; i++)
		out->size[i] = file->ncells_vec[i];
	
  start_4d[0] = block_no;
  start_4d[1] = 0;
  start_4d[2] = 0;
  start_4d[3] = 0;
  
  count_4d[0] = 1;
  count_4d[1] = file->ncells_vec[2];
  count_4d[2] = file->ncells_vec[1];
  count_4d[3] = file->ncells_vec[0];

  status = ncmpi_inq_varid(*((int*)file->handle), mangled_name, &varid);
  status = ncmpi_get_vara_double_all(*((int*)file->handle), varid, start_4d, count_4d, out->data);
  if (status < 0){
    printf("Error: ncmpi_get_vara_double_all, unknowns\n");
  }

  return out;
}

#endif
