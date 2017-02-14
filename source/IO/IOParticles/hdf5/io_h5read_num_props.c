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


void FTOC(io_h5read_num_props)(hid_t* file_identifier,
                                    int* num_part_props){

  hid_t dataset,dataspace, memspace;
  herr_t status;

  hsize_t dimens_2d[2], maxdimens_2d[2];

  dataset = H5Dopen(*file_identifier, "particle names");
  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
  
  *num_part_props = dimens_2d[0];
  
  
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  return;
}
