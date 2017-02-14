#define H5_USE_16_API

#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define NPART 16384
#define NPROP 17
#define POSX 9
#define VELX 14
#define FILE_NAME "pan2d_hdf5_chk_0004"
#define DATASET_NAME "tracer particles"
#define PLOT_NAME "particle.out"


/*
   Plots x-component of position against x-component of velocity.
   Use the parfile named test_paramesh_2d.par. Compile using
   a similar line to below:

   Looks at checkpoint file 4 which corresponds to zFinal = 0.
   Note, you will have to change zFinal in flash.par because the
   parfiles now go to zFinal=5 in coldstart parfile and then
   zFinal=0 in restart parfile.

   gcc validate_pancake.c \
   -I/opt/hdf5/1.8.4/openmpi-1.3.3_gnu-4.4.2/include \
   -I/opt/openmpi/1.3.3/gnu-4.4.2/include \
   -L/opt/hdf5/1.8.4/openmpi-1.3.3_gnu-4.4.2/lib -lhdf5 \
   -o validate_pancake

   After executing the program, you may visualize the
   data using gnuplot:
   plot "particle.out" using 1:2
*/

int main()
{
  hid_t fileID, dataset;
  double *pParticlesBuf;
  herr_t err;
  int i;
  FILE *file;

  fileID = H5Fopen(FILE_NAME, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(fileID >= 0);

  dataset = H5Dopen(fileID, DATASET_NAME);
  assert(dataset >=0);

  pParticlesBuf = malloc(NPART * NPROP * sizeof(*pParticlesBuf));

  err = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,
		H5P_DEFAULT, pParticlesBuf);
  assert(err >= 0);

  err = H5Dclose(dataset);
  assert(err >= 0);

  err = H5Fclose(fileID);
  assert(err >= 0);

  file = fopen(PLOT_NAME, "w");
  for (i=0; i<NPART; ++i) {
    fprintf(file, "%e %e\n", pParticlesBuf[POSX+(i*NPROP)],
	    pParticlesBuf[VELX+(i*NPROP)]);
  }
  fclose(file);
  free(pParticlesBuf);

  return 0;
}
