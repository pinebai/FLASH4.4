
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5write_blk_particle_info)(int* myPE, 
                    hid_t* file_identifier,
                    int blk_particle_info[][5],
                    int* local_blocks,
                    int* total_blocks,
                    int* global_offset)
{
  hid_t dataspace, dataset, memspace, dxfer_template;
  herr_t status;

  int rank;
  hsize_t dimens_2d[2];

  hsize_t start_2d[2];
  hsize_t stride_2d[2], count_2d[2];

  /* set the dimensions of the dataset */
  rank = 2;
  dimens_2d[0] = *total_blocks;
  dimens_2d[1] = 5;
  

  dataspace = H5Screate_simple(rank, dimens_2d, NULL);
  if(dataspace < 0) {
    Driver_abortFlashC("io_h5write_blk_particle_info: error can't create dataspace");
  }


  if ((*myPE == MASTER_PE) && (*global_offset != 0)) {
    dataset = H5Dopen(*file_identifier, "blk unique id"); 
    if(dataset < 0) {
      Driver_abortFlashC("Error: H5Dopen io_h5write_blk_unique_id\n");
    }
  }else {
    /* create the dataset */
    dataset = H5Dcreate(*file_identifier, "blk particle info", H5T_NATIVE_INT,
                  dataspace, H5P_DEFAULT);
  }    

  /* create the hyperslab -- this will differ on the different processors */
  start_2d[0] = (hsize_t) (*global_offset);
  start_2d[1] = 0;

  stride_2d[0] = 1;
  stride_2d[1] = 1;

  count_2d[0] = (hsize_t) (*local_blocks);
  count_2d[1] = 5;

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
                        stride_2d, count_2d, NULL);
  if(status < 0) {
    printf("%s\n","Error: Unable to select hyperslab for blk particle info dataspace");
    Driver_abortFlashC("Error: Unable to select hyperslab for blk particle info dataspace");
  }


  /* create the memory space */
  rank = 2;
  dimens_2d[0] = *local_blocks;
  dimens_2d[1] = 5;

  memspace = H5Screate_simple(rank, dimens_2d, NULL);

  /*set up for collective mode*/
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);
#ifdef IO_HDF5_PARALLEL
  if(HDF5_MODE == COLLECTIVE)
    H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
#endif

  /* write the data */
  status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                dxfer_template, blk_particle_info);

  if(status < 0) {
    printf("%s\n","Unable to write blk_particle_info");
    Driver_abortFlashC("Unable to write blk_particle_info");
  }
 

  H5Pclose(dxfer_template);
  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

}



