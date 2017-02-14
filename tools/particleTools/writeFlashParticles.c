#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "hdf5.h"
#include "flash_ptio.h"



int writeFlashParticles(int fid, double *localParticles, int localnp,
			int numPartProps, int offset){

  return writeFlashParticlesHDF5(fid, localParticles, localnp, 
				 numPartProps, offset);
  /* hid_t file_id = (hid_t) fid;
    
    hid_t dataspace, dataset, memspace;
    herr_t status;
    hsize_t dimens_2d[2], start_2d[2], count_2d[2], stride_2d[2];

    int rank, i;
    
    dataset = H5Dopen(file_id, "tracer particles");

    dataspace = H5Dget_space(dataset);

    start_2d[0] = (hsize_t) offset;
    start_2d[1] = (hsize_t) 0;

    stride_2d[0] = (hsize_t) 1;
    stride_2d[1] = (hsize_t) 1;

    count_2d[0] = (hsize_t) localnp;
    count_2d[1] = (hsize_t) numPartProps;

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
				 stride_2d, count_2d, NULL);
    
    rank = 2;
    dimens_2d[0] = (hsize_t) localnp;
    dimens_2d[1] = (hsize_t) numPartProps;


    memspace = H5Screate_simple(rank, dimens_2d, NULL);

    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
		      H5P_DEFAULT, localParticles);

    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
    
 
    return NORMAL_STATUS;*/
}
