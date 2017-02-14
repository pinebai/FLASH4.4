#include <stdlib.h>
#include <string.h>
#include "mfhdf.h"
#include "flash_reader.h"

/* Slurp: read a full dataset into target. Returns 1 on error, 0 on success */
int FR_slurp_HDF4(int32 handle, char *name, void *target){
  int32 dataset, index, rank;
  intn status;
  int32 dummyData_type, dummyNum_attrs; /* maybe NULL works? */  
  int32 dimsizes[MAX_VAR_DIMS], start[MAX_VAR_DIMS];
  int i;

  for(i=0; i<MAX_VAR_DIMS; i++)
    start[i]=(int32) 0;

  index = SDnametoindex(handle, name);
  if(index==FAIL){
    printf("Error reading \"%s\": dataset not found\n", name);
    return 1; 
  }

  dataset = SDselect(handle, index);
  status = SDgetinfo(dataset, NULL, &rank, dimsizes, 
		      &dummyData_type, &dummyNum_attrs);
  
  status = SDreaddata(dataset, start, NULL, dimsizes, target);
  if (status==FAIL) {
    printf("Error reading \"%s\"\n", name);
    return 1; 
  }

  SDendaccess(dataset);
  return 0;
}  

FR_File *FR_open_HDF4(char *filename){
  FR_File *out;
  int32 handle;
  char temp_varnames[FR_MAXVARS*FR_VAR_STRING_SIZE];
  int i;

  out = (FR_File *) malloc(sizeof(FR_File));
  out->handle = malloc(sizeof(int32));
  out->format = FR_HDF4;
  strcpy(out->filename, filename);
  
  /* Open file */
  handle = SDstart(filename, DFACC_READ);
  *((int32 *)out->handle) = handle;
  
  if(handle==FAIL)
    return NULL;

  /* Slurp header data */
  if(FR_slurp_HDF4(handle, "total blocks", &(out->nblocks)))
    return NULL;

  /* number of zones per block (dataset incorrectly named)*/
  if(FR_slurp_HDF4(handle, "number of blocks per zone", &(out->ncells_vec)))
    return NULL;

  /* use ncells_vec to get dim */
  if (out->ncells_vec[2] != 1)
    out->dim = 3;
  else if (out->ncells_vec[1] != 1)
    out->dim = 2;
  else
    out->dim = 1;
  
  out->ncells = out->ncells_vec[0] *
                out->ncells_vec[1] *
                out->ncells_vec[2] ;


  out->nodetype = (int *) malloc(sizeof(int)*out->nblocks);
  if(FR_slurp_HDF4(handle, "node type", out->nodetype))
    return NULL;

  out->lref = (int *) malloc(sizeof(int)*out->nblocks);
  if(FR_slurp_HDF4(handle, "refine level", out->lref))
    return NULL;

  out->coord = (double *) malloc(sizeof(double)*out->dim*out->nblocks);
  if(FR_slurp_HDF4(handle, "coordinates", out->coord))
    return NULL;

  out->size = (double *) malloc(sizeof(double)*out->dim*out->nblocks);
  if(FR_slurp_HDF4(handle, "block size", out->size))
    return NULL;

  out->bbox = (double *) malloc(sizeof(double)*out->dim*out->nblocks*2);
  if(FR_slurp_HDF4(handle, "bounding box", out->bbox))
    return NULL;

  /* now get variable names and number of variables */
  for(i=0; i<(FR_MAXVARS*FR_VAR_STRING_SIZE); i++)
    temp_varnames[i]='\0';

  if(FR_slurp_HDF4(handle, "unknown names", temp_varnames))
    return NULL;

  out->nvar = 0;
  for(i=0; i<FR_MAXVARS; i++){
    if(temp_varnames[i*FR_VAR_STRING_SIZE]=='\0')
      break;
    (out->nvar)++;
    strncpy(out->varnames[i], temp_varnames+FR_VAR_STRING_SIZE*i, 
	    FR_VAR_STRING_SIZE);
    out->varnames[i][FR_VAR_STRING_SIZE]='\0';
  }

  return out;
}

FR_Block *FR_GetBlock_HDF4(FR_File *file, char *var, int block_no){
  FR_Block *out;
  
  int32 dataset, index;
  int32 start[5], stride[5], edge[5]; /*5==1+MDIM+1 */
  int i, ivar;

  out = (FR_Block *) malloc(sizeof(FR_Block));
  out->data = (double *) malloc(sizeof(double)*file->ncells);

  /* this is fast enough */
  for(ivar=0; ivar<file->nvar; ivar++)
    if(strcmp(var, file->varnames[ivar])==0)
      break;

  /* unk is unk[block_no, nzb, nyb, nxb, ivar] here */
  start[0]=block_no;
  for(i=1; i<4; i++)
    start[i]=0;
  start[4]=ivar;
  
  edge[0]=1;
  edge[4]=1;
  for(i=1; i<4; i++)
   edge[i]=file->ncells_vec[3-i];

  for(i=0; i<5; i++)
    stride[i]=1;

  index = SDnametoindex(*((int32 *)file->handle), "unknowns");
  dataset = SDselect(*((int32 *)file->handle), index);
  SDreaddata(dataset, start, stride, edge, out->data);
  SDendaccess(dataset);

  return out;
}
