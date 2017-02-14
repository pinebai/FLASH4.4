#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
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
void FTOC(io_h5read_facevars)(hid_t* file_identifier,
                   int* globalNXB,  /* num of zones in x over whole domain*/
                   int* globalNYB,  /* num of zones in y over whole domain*/
                   int* globalNZB,  /* num of zones in z over whole domain*/
                   int *nxbOffset,  /* num of zones in x right of 0 to read*/
		   int *nybOffset,  /* num of zones in y right of 0 to read*/
		   int *nzbOffset,  /* num of zones in z right of 0 to read*/
		   int *nxb,  /*num of zones in x to be read*/
		   int *nyb,  /*num of zones in y to be read*/
		   int *nzb,  /*num of zones in z to be read*/
		   double* unknowns,     /* [mblk][NZB][NYB][NXB][var] */
                   char record_label[5]) /* add char--null termination */
{
  hid_t dataspace, dataset, memspace;
  herr_t status;

  int rank;
  hsize_t dimens_5d[5];

  hsize_t start_4d[4], start_5d[5];
  hsize_t stride_4d[4], stride_5d[5], count_4d[4], count_5d[5];

  char record_label_new[5];

  int ierr;
  int total_blocks = 1; /* doesn't mean all that much for nofbs */

  /* the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination */
  strncpy(record_label_new, record_label,4);

  *(record_label_new + 4) = '\0';

 /* if(nxbOffset > 0)
    nxbOffset++;
  if(nybOffset > 0)
    nybOffset++;
  if(nzbOffset > 0)
    nzbOffset++;
*/

  /* open the dataset */
  dataset = H5Dopen(*file_identifier, record_label_new);

  if (dataset < 0) {
    printf("couldn't find variable '%s' in the file, so skipping it\n", record_label_new);
  } else {
    dataspace = H5Dget_space(dataset);

    /* create the hyperslab -- this will differ on the different processors */
    start_4d[0] = 0; /*(hsize_t) (*global_offset);*/
    start_4d[1] = *nzbOffset;
    start_4d[2] = *nybOffset;
    start_4d[3] = *nxbOffset;

    stride_4d[0] = total_blocks;
    stride_4d[1] = 1;
    stride_4d[2] = 1;  
    stride_4d[3] = 1;

    count_4d[0] = total_blocks;
    count_4d[1] = *nzb;
    count_4d[2] = *nyb;
    count_4d[3] = *nxb;

    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
                         stride_4d, count_4d, NULL);
  
    if (status < 0){
      printf("Error: Unable to select hyperslab for unknowns dataspace\n");
      Driver_abortFlashC("Error: Unable to select hyperslab for unknowns dataspace\n");
    }

    /* create the memory space */
    rank = 5;
    dimens_5d[0] = total_blocks;
    dimens_5d[1] = *nzb;
    dimens_5d[2] = *nyb;
    dimens_5d[3] = *nxb;
    dimens_5d[4] = 1;

    memspace = H5Screate_simple(rank, dimens_5d, NULL);

    /*
    start_5d[0] = 0;
    start_5d[1] = (*nguard)*k3d;
    start_5d[2] = (*nguard)*k2d;
    start_5d[3] = (*nguard);
    start_5d[4] = *index-1;  

    stride_5d[0] = 1;
    stride_5d[1] = 1;
    stride_5d[2] = 1;
    stride_5d[3] = 1;
    stride_5d[4] = 1;

    count_5d[0] = *local_blocks;
    count_5d[1] = *nzb;
    count_5d[2] = *nyb;
    count_5d[3] = *nxb;
    count_5d[4] = 1;
   
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                         start_5d, stride_5d, count_5d, NULL);
    

    if (ierr < 0){
      printf("Error: Unable to select hyperslab for coordinates memspace\n");
      Driver_abortFlash("Error: Unable to select hyperslab for coordinates memspace\n");
    }*/


    /* read the data */
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                 H5P_DEFAULT, unknowns);
               
    if (status < 0){
      printf("Error: Unable to read unknowns\n %d\n",status );
      Driver_abortFlashC("Error: Unable to read unknowns\n");
    }

    H5Sclose(memspace); 
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
}
