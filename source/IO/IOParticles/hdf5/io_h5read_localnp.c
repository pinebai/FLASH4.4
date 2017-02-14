#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Flash.h"
#include "constants.h"




int Driver_abortFlashC(char* message);


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/* This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from the varnames or species
   database */
void FTOC(io_h5read_localnp)(hid_t* file_identifier,
                       int* localnp,     /* [mblk][NZB][NYB][NXB][var] */
                       int* local_blocks,
                       int* total_blocks,
                       int* global_offset)
{
  hid_t dataspace, dataset, memspace;
  herr_t status;
  
  int rank;
  
  hsize_t start_1d, stride_1d, count_1d, dimens_1d;
  
  int ierr;
  

  
  /* open the dataset */
  dataset = H5Dopen(*file_identifier, "localnp");
  if(dataset < 0){
    Driver_abortFlashC("Error, unable to read dataset localnp");
  }


  dataspace = H5Dget_space(dataset);
  if(dataspace < 0){
    Driver_abortFlashC("Error, unable to read dataset localnp");
  }


  /* create the hyperslab -- this will differ on the different processors */
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_blocks);
  

  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
                         &stride_1d, &count_1d, NULL);
  
  if (status < 0){
    printf("Error: Unable to select hyperslab for localnp dataspace\n");
    Driver_abortFlashC("Error: Unable to select hyperslab for localnp dataspace\n");
  }
  

  /* create the memory space */
  rank = 1;
  dimens_1d = *local_blocks;
  
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  if(memspace < 0){
    Driver_abortFlashC("Error: Unable to create memspace\n");
  }


  /* read the data */
  status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, 
               H5P_DEFAULT, localnp);
  
  if (status < 0){
    printf("Error: Unable to read localnp\n");
    Driver_abortFlashC("Error: Unable to read localnp\n");
  }

  
  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);


}

























































