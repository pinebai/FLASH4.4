#include <stdio.h>
#include <stdlib.h>
#include "flash_reader.h"

/* test_reader(filename)

   Tries out the FR_ routines to see if they work. Would be better if it came
   with a sample hdf file, we could test stuff like number precision.
*/

void test_reader(char *filename){
  FR_File *file;
  FR_Block *block;
  int i,j,k,d;
  FR_ParticlesAllProps *particlesAllProps;
  int particlesRead, particlesGotten;

  printf("opening: %s\n", filename);
  file = FR_open(filename);

  if(file->format==FR_HDF5)
    printf("HDF5\n");
  else if(file->format==FR_HDF4)
    printf("HDF4\n");
  else if (file->format==FR_NCDF)
    printf("NCDF\n");
  else
    printf("Mystery format!\n");

  printf("filename: %s\nnblocks: %d\n", file->filename, file->nblocks);
 
  printf("block shape: ");
  for(i=0;i<3;i++)
    printf("%d ", file->ncells_vec[i]);

  printf("\nncells: %d\n", file->ncells);

  printf("dim: %d\n", file->dim);
  printf("nodetype:\n");
  for(i=0;i<file->nblocks;i++)
    printf("%d ", file->nodetype[i]);
  printf("\n");
#if 1
  printf("coords:\n");
  for(i=0;i<file->nblocks;i++){
    for (d=0;d<file->dim;d++)
      printf("%g ", file->coord[3*i+d]);
    printf("\n");
    }

  printf("\nblock size:\n"); 
  for(i=0;i<file->nblocks;i++){
    for (d=0;d<file->dim;d++)
      printf("%g ", file->size[3*i+d]);
    printf("\n");
    }
#endif  
  printf("nvar: %d\n",file->nvar);

  for(i=0;i<file->nvar;i++)
    printf("%s\n", file->varnames[i]);


  printf("Trying to read in blocks...\n");
  for(i=0; i<file->nvar; i++){
    printf("%s: ", file->varnames[i]);
    for(j=0; j<file->nblocks; j++){
      printf("%d ", j);
      block = FR_GetBlock(file, file->varnames[i], j);
      for(k=0; k<file->ncells; k++)
	block->data[k] = 0.; /* Checking for bad memory */
      FR_DeleteBlock(block);
    }
    printf("\n");
  }
  printf("done reading blocks\n");
#if 1
  printf("dumping block %d, %s\n", 1, file->varnames[0]);
  block = FR_GetBlock(file, file->varnames[0], 1);
  for(k=0; k<file->ncells; k++)
   printf("%g ", block->data[k]);
  FR_DeleteBlock(block);
  printf("\n");
#endif

  particlesRead = 0;

  particlesRead = 0;
  printf("Trying to read in particlesAllProps...\n");
  while (particlesRead < file->totalparticles) {
    particlesAllProps = FR_GetParticlesAllProps(file, particlesRead, 100, &particlesGotten);
    if (particlesRead == 0) {
      printf("dumping first set of particles\n");
      for (j = 0; j < particlesAllProps->numParticles; j++) {
	printf("particle: %d\n", j);
	printf(" realprops:\n");
	for (i = 0; i < particlesAllProps->numRealProps; i++) {
	  printf("  %24s: %f\n", particlesAllProps->realPropsNames[i], particlesAllProps->realProps[j*(particlesAllProps->numRealProps) + i]);
	}
	printf(" intprops:\n");
	for (i = 0; i < particlesAllProps->numIntProps; i++) {
	  printf("  %24s: %d\n", particlesAllProps->intPropsNames[i], particlesAllProps->intProps[j*(particlesAllProps->numIntProps) + i]);
	}
      }
    }

    particlesRead += particlesGotten;
    printf("Read %d; total %d\n", particlesGotten, particlesRead);
    FR_DeleteParticlesAllProps(particlesAllProps);
  }
  printf("done reading particles\n");

  FR_close(file);
}
