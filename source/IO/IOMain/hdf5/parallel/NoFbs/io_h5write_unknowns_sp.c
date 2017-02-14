#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Flash.h"
#include "constants.h"
#include "io_flash.h"
#include "io_h5_attribute.h"

int Driver_abortFlashC(char* message);

/* 
   This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from the varnames or species
   database 
   
   The dimensions of the unknowns array (nvar, nxb, nyb, nzb, maxblocks)
   are passed through as arguments.  The dataspace (what gets written
   to disk) and the memory space (how the unknowns array is accessed in
   local memory) are defined based on these passed values.  This allows
   use to pass an unk array that has all the guardcells + interior cells
   (as in the checkpointing), or a dummy unk array that has just the 
   interior cells (in this case, nguard would be passed as 0).
*/


void FTOC(io_h5write_unknowns_sp)(
  int* myPE,
  hid_t* file_identifier,
  int* globalNXB,       /* num of zones in x dir across whole domain */
  int* globalNYB,       /* num of zones in y dir across whole domain */
  int* globalNZB,       /* num of zones in z dir across whole domain */
  int* nxbOffset,       /* num of zones in x dir to right of 0 where we'll start writing */
  int* nybOffset,       /* num of zones in y dir to right of 0 where we'll start writing */
  int* nzbOffset,       /* num of zones in z dir to right of 0 where we'll start writing */
  int* nxb,             /* num of zones in x dir to be written on this call */
  int* nyb,             /* num of zones in y dir to be written on this call */
  int* nzb,             /* num of zones in z dir to be written on this call */
  float* unknowns,     /* [1][NZB][NYB][NXB] */
  float* varMin,       /* minimum value of unknown */
  float* varMax,       /* maximum value of unknown */
  char record_label[5],/* add char-null termination */
  int* dowrite         /* T/F (1/0) bool to send data */
) {

  hid_t dataspace, dataset, memspace, dxfer_template, dataset_plist;

  herr_t status, ierr;

  int rank;
  hsize_t dimens_4d[4], dimens_5d[5];

  hsize_t start_4d[4];
  hsize_t stride_4d[4], count_4d[4];

#ifdef CHUNK
  hsize_t dimens_chunk[4];
#endif

  char record_label_new[5];

  int total_blocks = 1; /* not necessary for nofbs, but keeps consistency for fidlr routines */

  const char minLabel[] = "minimum";
  const char maxLabel[] = "maximum";
  const int attType = IO_FLASH_FLOAT;
  const int dims = 1;
  const int diskSize[] = {1};

  /*printf("varMin = %e\nvarMax = %e\n", *varMin, *varMax);*/
  /*  printf("myPE: %d\nglobalNXB: %d\nglobalNYB: %d\nglobalNZB: %d\nnxbOffset: %d\nnybOffset: %d\nnzbOffset%d\nnxb: %d\nnyb: %d\nnzb: %d\n",
   *myPE, *globalNXB, *globalNYB, *globalNZB, *nxbOffset, *nybOffset, *nzbOffset, *nxb, *nyb, *nzb);*/


  /* 
     the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination 
  */
  strncpy(record_label_new, record_label,4);
  *(record_label_new + 4) = '\0';

  /* set the dimensions of the dataset */
  rank = 4;
  dimens_4d[0] = total_blocks;
  dimens_4d[1] = *globalNZB;
  dimens_4d[2] = *globalNYB;
  dimens_4d[3] = *globalNXB;

  dataspace = H5Screate_simple(rank, dimens_4d, NULL);
  if(dataspace < 0){
    printf("io_h5write_unknowns: dataspace error");
    Driver_abortFlashC("io_h5write_unknowns: dataspace error");
  }

  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);
  if(dataset_plist < 0){
    printf("io_h5write_unknowns: dataset_plist error");
    Driver_abortFlashC("io_h5write_unknowns: dataset_plist error");
  }

    
  /* create a parallel hdf5 file */
  dataset = H5Dcreate(*file_identifier, record_label_new,
	              H5T_NATIVE_FLOAT, dataspace, dataset_plist); 
  if(dataset < 0) {
    Driver_abortFlashC("dataset Error: H5Dcreate io_h5write_unk\n");
  }    


  start_4d[0] = 0;
  start_4d[1] = *nzbOffset; /*number of z zones to the 'left' of myPE -- do this by copying functionality in Particles_getOffset */
  start_4d[2] = *nybOffset;
  start_4d[3] = *nxbOffset;

  stride_4d[0] = total_blocks;
  stride_4d[1] = 1;
  stride_4d[2] = 1;  
  stride_4d[3] = 1;

  count_4d[0] = total_blocks;
  count_4d[1] = *nzb;  /*now really write the local nxb, nyb and nzb */
  count_4d[2] = *nyb;  /*remember still 1 block per proc, it's just that the blocks can have various, nxb, nyb, nzb */
  count_4d[3] = *nxb;

  ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
			     stride_4d, count_4d, NULL);

  if(ierr < 0){
     printf("%s\n", "Error: unable to select hyperslab for unknowns dataspace");
     Driver_abortFlashC("Error: unable to select hyperslab for unknowns dataspace");
  }


  rank = 5;
  dimens_5d[0] = total_blocks;
  dimens_5d[1] = *nzb;
  dimens_5d[2] = *nyb;
  dimens_5d[3] = *nxb;
  dimens_5d[4] = 1;


  memspace = H5Screate_simple(rank, dimens_5d, NULL);
  if(memspace < 0){
    printf("io_h5write_unknowns: memspace error");
    Driver_abortFlashC("io_h5write_unknowns: memspace error");
  }

  /* obtain a copy of the file transfer property list */ 
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);


  /* default for now */
  /* default is 'independent' mode */


  /*ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);*/
  /*ierr = H5Pset_preserve(dxfer_template, 0u);*/

 


#ifdef CHUNK
  /* set the layout to chunked */

  ierr = H5Pset_layout(dataset_plist, H5D_CHUNKED);
  
  /* create a chunk a containing 10 blocks worth of data */
  dimens_chunk[0] = 10;
  dimens_chunk[1] = *nzb;
  dimens_chunk[2] = *nyb;
  dimens_chunk[3] = *nxb;
   
  ierr = H5Pset_chunk(dataset_plist, 4, dimens_chunk);
#endif


  /*!!DEV -kda I've removed the max and min parts for now */

  /* write the data */
  if(!*dowrite) {
    ierr = H5Sselect_none(dataspace);
    if(ierr < 0) Driver_abortFlashC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
    ierr = H5Sselect_none(memspace);
    if(ierr < 0) Driver_abortFlashC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
  }
  status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, dxfer_template, unknowns);
  if(status < 0) {
    printf("io_h5write_unknowns: H5Dwrite error");
    Driver_abortFlashC("io_h5write_unknowns: H5Dwrite error");
  }
#ifdef DEBUG_IO
  printf("UNKNOWNS: wrote unknowns, status = %d\n", (int) status);
#endif
  
  H5Pclose(dxfer_template);

  H5Pclose(dataset_plist);

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);

  /* Write attributes for unknowns: varMin and varMax.
     Unlike fixed block size mode, we only have to worry about this in
     parallel for now.  Note: The io_h5_attribute_XXX functions open
     the dataset "record_label_new" before creating / writing
     attributes, so "record_label_new" should be closed before the
     call. */
  io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			 diskSize, record_label_new, minLabel);
  io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			record_label_new, minLabel, varMin);

  io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			 diskSize, record_label_new, maxLabel);
  io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			record_label_new, maxLabel, varMax);
}
