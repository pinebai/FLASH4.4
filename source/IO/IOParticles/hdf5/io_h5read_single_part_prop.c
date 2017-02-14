#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"

int Driver_abortFlashC(char* message);


void FTOC(io_h5read_single_part_prop)(hid_t* file_identifier,
                                    double *particles,
                                    int *localnp,
                                    int *part_prop_id,
                                    int *file_part_prop,
                                    int *file_npart_props,
                                    int *current_npart_props,
                                    int *particle_offset)
{

  hid_t    dataspace, dataset, memspace, xfer_plist;
  herr_t   status, herr;

  int      rank;
  hsize_t  dimens_2d[2], start_2d[2], count_2d[2], stride_2d[2];

  if(*localnp <= 0){
    /*we shouldn't have to do anything here*/
    return;
  }

  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef IO_HDF5_PARALLEL
  if (HDF5_MODE == COLLECTIVE) {
    herr = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(herr >= 0);
  }
#endif


  dataset = H5Dopen(*file_identifier, "tracer particles");
  if(dataset < 0){
    Driver_abortFlashC("[IO ERROR]: Could not open tracer particles for reading!\n");
  }

  dataspace = H5Dget_space(dataset);

  
  /*for dataspace - this is the file's end of things*/
  start_2d[0] = (hsize_t)(*particle_offset);
  start_2d[1] = (hsize_t)(*file_part_prop - 1);

  stride_2d[0] = 1;
  stride_2d[1] = (hsize_t)(*file_npart_props);
  
  count_2d[0] = (hsize_t)(*localnp);
  count_2d[1] = 1;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                               stride_2d, count_2d, NULL);

  if (status < 0){
    printf("Error: Unable to select hyperslab for particles dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for particles dataspace\n");
  }


  /*for memspace - this is the current simulation's end of things*/
  rank = 2;
  dimens_2d[0] = (hsize_t)(*localnp);
  dimens_2d[1] = (hsize_t)(*current_npart_props);

  start_2d[0] = (hsize_t)0;
  start_2d[1] = (hsize_t)(*part_prop_id - 1);

  stride_2d[0] = 1;
  stride_2d[1] = (hsize_t)(*current_npart_props);  

  count_2d[0] = (hsize_t)(*localnp);
  count_2d[1] = 1;

  memspace = H5Screate_simple(rank, dimens_2d, NULL);

  status = H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start_2d,
                               stride_2d,count_2d, NULL);
  if (status < 0){
    printf("Error: Unable to select hyperslab for particles memspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for particles memspace\n");
  }
  
  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                   xfer_plist, particles);
  if (status < 0){
    H5Eprint(stdout);
    printf("[io_h5read_single_part_prop] Error: Unable to read particles from dataset\n");
    Driver_abortFlashC("[io_h5read_single_part_prop] Error: Unable to read particles from dataset\n");
  }

  H5Pclose(xfer_plist);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return;
  
}
