#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5write_generic_real_arr)(int* myPE, 
				       hid_t* file_identifier,
				       double* generic_arr,
				       int* local_size,
				       int* total_size,
				       int* global_offset,
				       char* dataset_name, 
				       int* name_len)
{
  hid_t dataspace, dataset, memspace, dataset_plist;
  herr_t status;

  int rank;
  hsize_t dimens_1d;

  hsize_t start_1d;
  hsize_t stride_1d, count_1d;


  char* dataset_name_new;

  dataset_name_new = (char *) malloc((*name_len) + 1 * sizeof(char)); 

  /* copy the dataset_name into a c string dataset_name_new with 
     its exact length with the \0 termination */

  strncpy(dataset_name_new, dataset_name, *name_len);
  *(dataset_name_new + *name_len) = '\0';



  /* set the dimensions of the dataset */
  rank = 1;
  dimens_1d = *total_size;
 
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_size);

 
  dataspace = H5Screate_simple(rank, &dimens_1d, NULL);
  if(dataspace < 0) {
     Driver_abortFlashC("Error: H5Screate_simple io_h5write_generic_arr\n");
  }



  /*DEV: this line was only in the serial version */
  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);
  if(dataset_plist < 0) {
    Driver_abortFlashC("Error: dataset_plist io_h5write_generic_arr\n");
  }
  

  /*This part is necessary when serial IO is used
    The master proc writes all the data but only creates
    the dataset once.  In the parallel IO case each 
    proc calls H5Dcreate*/
  if ((*myPE == MASTER_PE) && (*global_offset != 0)) {
     dataset = H5Dopen(*file_identifier, dataset_name_new); 
    if(dataset < 0) {
       Driver_abortFlashC("Error: H5Dopen io_h5write_generic_arr\n");
    }
   
  }else {

     /* create the dataset */
    dataset = H5Dcreate(*file_identifier, dataset_name_new, H5T_NATIVE_DOUBLE,
                  dataspace, dataset_plist);
    if(dataset < 0) {
       Driver_abortFlashC("Error: H5Dcreate io_h5write_generic_arr\n");
    }

    /*dataset_plist was H5P_DEFAULT*/
  }



  if(*local_size > 0) {

    /* create the hyperslab -- this will differ on the different processors */
    
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
				 &stride_1d, &count_1d, NULL);
    
    if(status < 0) {
      Driver_abortFlashC("Error: H5Sselect_hyperslab io_h5write_generic_arr\n");
    }
    
    
    /* create the memory space */
    dimens_1d = *local_size;
    memspace = H5Screate_simple(rank, &dimens_1d, NULL);
    if(memspace < 0) {
      Driver_abortFlashC("Error: H5Screate_simple mem io_h5write_generic_arr\n");
    }
  

    
    
    
    /* write the data */
    if(*local_size == *total_size){
      if(*myPE == MASTER_PE){
	status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
			  H5P_DEFAULT, generic_arr);
    
      }
    }else{
      status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
			H5P_DEFAULT, generic_arr);
    }
 
    if(status < 0) {
      Driver_abortFlashC("Error: H5Dwrite io_h5write_generic_real_arr\n");
    }
    

    H5Sclose(memspace); 
    H5Pclose(dataset_plist);
  }
    


  H5Sclose(dataspace);
  H5Dclose(dataset);

  free(dataset_name_new);
  
}


