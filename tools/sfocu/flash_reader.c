#include <stdlib.h>
#include "flash_reader.h"
#include "options.h"

#ifndef NO_HDF4
#include "mfhdf.h"
extern FR_File *FR_open_HDF4(char *filename);
extern FR_Block *FR_GetBlock_HDF4(FR_File *file, char *var, int block_no);
#endif

#ifndef NO_HDF5
#include "hdf5.h"
extern FR_File *FR_open_HDF5(char *filename, options_t opts);
extern FR_File *FR_open_HDF5_Pmesh(hid_t handle, FR_File* out, options_t opts);
extern FR_File *FR_open_HDF5_Chombo(hid_t handle, FR_File* out, options_t opts);
extern FR_Block *FR_GetBlock_HDF5(FR_File *file, char *var, int block_no);
extern FR_Block *FR_GetBlock_HDF5_Chombo(FR_File *file, int var, int block_no);
extern FR_ParticlesAllProps *FR_GetParticlesAllProps_HDF5(FR_File *file, int startingAt, int atMost, int *numberGotten);
#endif

#ifndef NO_NCDF
#include "pnetcdf.h"
extern FR_File *FR_open_NCDF(char *filename, options_t opts);
extern FR_Block *FR_GetBlock_NCDF(FR_File *file, char *var, int block_no);
extern FR_ParticlesAllProps *FR_GetParticlesAllProps_NCDF(FR_File *file, int startingAt, int atMost, int *numberGotten);
#endif

FR_Block *FR_GetBlock(FR_File *file, char *var, int block_no){
  switch (file->format){

#ifndef NO_HDF4
  case FR_HDF4:
    return FR_GetBlock_HDF4(file, var, block_no);
    break;
#endif

#ifndef NO_HDF5
  case FR_HDF5_PMESH:
    return FR_GetBlock_HDF5(file, var, block_no);
    break;
#endif

#ifndef NO_NCDF
  case FR_NCDF:
    return FR_GetBlock_NCDF(file, var, block_no);
    break;
#endif

  default:
    return NULL;
  }
}

FR_ParticlesAllProps *FR_GetParticlesAllProps(FR_File *file, int startingAt, 
			     int atMost, int *numberGotten) {

  switch (file->format){

#ifndef NO_HDF4
  case FR_HDF4:
    /* return FR_GetParticles_HDF4(file, var, block_no); */
    break;
#endif

#ifndef NO_HDF5
  case FR_HDF5_PMESH:
    return FR_GetParticlesAllProps_HDF5(file,startingAt, atMost, numberGotten);
    break;
#endif

#ifndef NO_NCDF
  case FR_NCDF:
    return FR_GetParticlesAllProps_NCDF(file,startingAt, atMost, numberGotten);
    break;
#endif

  default:
    return NULL;
  }
}

void FR_DeleteBlock(FR_Block *block){
  free(block->data);
  free(block);
}

void FR_DeleteParticlesAllProps(FR_ParticlesAllProps *particles){
  int i;
  for (i = 0; i < particles->numRealProps; i++) {
    free(particles->realPropsNames[i]);
  }

  for (i = 0; i < particles->numIntProps; i++) {
      free(particles->intPropsNames[i]);
  }
    
  free(particles->intPropsNames);
  free(particles->realPropsNames);
  free(particles->realProps);
  free(particles->intProps);
  free(particles);
}

FR_File *FR_open(char *filename, options_t opts){
  int c, format;
  int isNCDF, isHDF5;
  FILE *fp;
  int ch[8];

  /* From H5Fpkg.h */
  /* const char *H5F_SIGNATURE = "\211HDF\r\n\032\n"; */
  const int H5F_SIGNATURE[8] = {137, 72, 68, 70, 13, 10, 26, 10};
  const int H5F_SIGNATURE_LEN = 8;

  const char *NCDF_SIGNATURE = "CDF";
  const int NCDF_SIGNATURE_LEN = 3;

  if((fp = fopen(filename, "r"))){ 
    for(c=0; c<H5F_SIGNATURE_LEN; c++) {
      ch[c] = fgetc(fp);
    }
  
    isHDF5 = 1;
    isNCDF = 1;
    for (c=0; c<H5F_SIGNATURE_LEN; c++) {
      if (ch[c] != H5F_SIGNATURE[c]) {
	isHDF5 = 0;
	break;
      }
    }
  
    for (c=0; c<NCDF_SIGNATURE_LEN; c++) {
      if (ch[c] != NCDF_SIGNATURE[c]) {
	isNCDF = 0;
	break;
      }
    }
    
    if (isHDF5 == 1) {
      format = FR_HDF5;
    } else if (isNCDF == 1) {
      format = FR_NCDF;
    } else {
      format = FR_HDF4;
    }
    fclose(fp);
  }

  switch (format){
#ifndef NO_HDF4
  case FR_HDF4:
    return FR_open_HDF4(filename);
    break;
#endif

#ifndef NO_HDF5
  case FR_HDF5:
    return FR_open_HDF5(filename, opts);
    break;
#endif

#ifndef NO_NCDF
  case FR_NCDF:
    return FR_open_NCDF(filename, opts);
    break;
#endif

  default:
    return NULL;
  }
}

void FR_close(FR_File *file){
  free(file->nodetype);
  free(file->coord);
  free(file->size);
  free(file->bbox);
  free(file->lref);
  switch (file->format){

#ifndef NO_HDF4
  case FR_HDF4:
    SDend(*((int32 *)file->handle)); 
    break;
#endif

#ifndef NO_HDF5
  case FR_HDF5:
    H5Fclose(*((hid_t *)file->handle));
    break;
#endif

#ifndef NO_NCDF
  case FR_NCDF:
    ncmpi_close(*((int *)file->handle));
    break;
#endif

  default:
    break;
  }
  free(file->handle);
  free(file);
}
