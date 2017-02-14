#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include <mpi.h>
#include "Flash.h"
#include "constants.h"
#include "io_flash.h"
#include "io_h5_attribute.h"


int Driver_abortFlashC(char* message);


/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */



/* 
   This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from an unk variable in single precision
   
   
   The dimensions of the unknowns array (nxb, nyb, nzb, maxblocks)
   are passed through as arguments.  The dataspace (what gets written
   to disk) and the memory space (how the unknowns array is accessed in
   local memory) are defined based on these passed values. 
   We also pass in the maximum and minimum values of the dataset 

*/


void FTOC(io_h5write_unknowns_sp)(int* myPE, 
                          hid_t* file_identifier,
                          int* nxb,            /* # of zones to store in x */
                          int* nyb,            /* # of zones to store in y */
                          int* nzb,            /* # of zones to store in z */
                          float* varMin,
                          float* varMax,
                          float* unknowns,   /* [mblk][NZB][NYB][NXB][nvar] */
                          char record_label[5],/* add char-null termination */
                          int* local_blocks,  
                          int* total_blocks,
                          int* global_offset,
                          int* dowrite)
{


  hid_t dataspace, dataset, memspace, dxfer_template, dataset_plist;

  herr_t status;

  int rank;
  hsize_t dimens_4d[4], dimens_5d[5];

  hsize_t start_4d[4];
  hsize_t stride_4d[4], count_4d[4];

#if defined (CHUNK) || defined (IBM) || defined (TFLOPS) || defined (CHIBA) || defined (SGI)
  int ierr;
#endif

#ifdef CHUNK
  hsize_t dimens_chunk[4];
#endif

  char record_label_new[5];
  const char minLabel[] = "minimum";
  const char maxLabel[] = "maximum";
  const int attType = IO_FLASH_FLOAT;
  const int dims = 1;
  const int diskSize[] = {1};
  int parallelIO;
#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif


  /* 
     the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination 
  */
  strncpy(record_label_new, record_label,4);
  *(record_label_new + 4) = '\0';

  /* set the dimensions of the dataset */
  rank = 4;
  dimens_4d[0] = *total_blocks;
  dimens_4d[1] = *nzb;
  dimens_4d[2] = *nyb;
  dimens_4d[3] = *nxb;
  
  dataspace = H5Screate_simple(rank, dimens_4d, NULL);

#ifdef DEBUG_IO
  printf("UNKNOWNS: dataspace = %d\n", (int) dataspace);
#endif

  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);

#ifdef CHUNK
  /* set the layout to chunked */

  ierr = H5Pset_layout(dataset_plist, H5D_CHUNKED);

  if (ierr < 0){
    printf("Error: Unable to set layout for dataset_plist\n");
    Driver_abortFlashC("Error: Unable to set layout for dataset_plist\n");
  }

  
  /* create a chunk containing 10 blocks worth of data */
  dimens_chunk[0] = 10;
  dimens_chunk[1] = *nzb;
  dimens_chunk[2] = *nyb;
  dimens_chunk[3] = *nxb;
   
  ierr = H5Pset_chunk(dataset_plist, 4, dimens_chunk);

  if (ierr < 0){
    printf("Error: Unable to set chunk for dataset_plist\n");
    Driver_abortFlashC("Error: Unable to set chunk for dataset_plist\n");
  }

#endif


  if ((*myPE == MASTER_PE) && (*global_offset != 0)) {
    dataset = H5Dopen(*file_identifier, record_label_new); 
    if(dataset < 0) {
      Driver_abortFlashC("Error: H5Dopen io_h5write_unknowns_sp\n");
    }
  }else {
    
    dataset = H5Dcreate(*file_identifier, record_label_new,
                  H5T_NATIVE_FLOAT, dataspace, dataset_plist); 
    if(dataset < 0) {
      Driver_abortFlashC("Error: H5Dcreate io_h5write_unknowns_sp\n");
    }    
  }


  if (*local_blocks > 0) {
    /* create the hyperslab -- this will differ on the different processors */
    start_4d[0] = (hsize_t) (*global_offset);
    start_4d[1] = 0;
    start_4d[2] = 0;
    start_4d[3] = 0;
    
    stride_4d[0] = 1;
    stride_4d[1] = 1;
    stride_4d[2] = 1;  
    stride_4d[3] = 1;
    
    count_4d[0] = (hsize_t) (*local_blocks);
    count_4d[1] = *nzb;
    count_4d[2] = *nyb;
    count_4d[3] = *nxb;
    
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
                         stride_4d, count_4d, NULL);
    
    if (status < 0){
      printf("Error: Unable to select hyperslab for unknowns dataspace\n");
      Driver_abortFlashC("Error: Unable to select hyperslab for unknown dataspace\n");
    }


    /* create the memory space -- we can get away with a simple memory space
       for the unknowns, since we are passing a contiguous block of memory
       now, and the block count is the last index in FORTRAN */
    rank = 5;
    dimens_5d[0] = *local_blocks;
    dimens_5d[1] = *nzb;
    dimens_5d[2] = *nyb;
    dimens_5d[3] = *nxb;
    dimens_5d[4] = 1;
     
    memspace = H5Screate_simple(rank, dimens_5d, NULL);

    /* setting the transfer template */
    dxfer_template = H5Pcreate(H5P_DATASET_XFER);



    /* ---------------------------------------------------------------------
       platform dependent code goes here -- the access template must be
       tuned for a particular filesystem blocksize.  some of these 
       numbers are guesses / experiments, others come from the file system
       documentation.
       
       The sieve_buf_size should be equal a multiple of the disk block size
       ---------------------------------------------------------------------- */

    /*
#ifdef IBM
    ierr = H5Pset_hyper_cache(dxfer_template, 1u, 0u);
    ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
    ierr = H5Pset_btree_ratios(dxfer_template, 0.334, 0.333, 0.333); 

     for best performance, set this to 0 
    ierr = H5Pset_preserve(dxfer_template, 0u);  
#endif
    
#ifdef TFLOPS
    ierr = H5Pset_hyper_cache(dxfer_template, 1u, 0u);
    ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
    ierr = H5Pset_preserve(dxfer_template, 0u);  
#endif

#ifdef CHIBA
    ierr = H5Pset_hyper_cache(dxfer_template, 1u, 0u);
    ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
    ierr = H5Pset_preserve(dxfer_template, 0u);  
#endif

#ifdef SGI
    ierr = H5Pset_hyper_cache(dxfer_template, 1u, 0u);
    ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);
    ierr = H5Pset_preserve(dxfer_template, 0u);  
#endif

    */

    /* ----------------------------------------------------------------------
       end of platform dependent data-transfer property list settings
       ---------------------------------------------------------------------- */
    /*not platform dependent, just a switchable mode*/
#ifdef IO_HDF5_PARALLEL
    if( HDF5_MODE == COLLECTIVE)
      H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
#endif
    if(!*dowrite) {
      status = H5Sselect_none(dataspace);
      if(status < 0) Driver_abortFlashC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
      status = H5Sselect_none(memspace);
      if(status < 0) Driver_abortFlashC("[" FILE_AT_LINE_C "] H5Sselect_none invalid return.");
    }
    /* write the data */
    status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, 
                  dxfer_template, unknowns);
    
    if (status < 0){
      printf("Error: Unable to write unknowns\n");
      Driver_abortFlashC("Error: Unable to write unknowns\n");
    }

#ifdef DEBUG_IO
    printf("UNKNOWNS: wrote unknowns, status = %d\n", (int) status);
#endif

    H5Pclose(dxfer_template);
    H5Sclose(memspace); 
  }  
  H5Pclose(dataset_plist);
  H5Sclose(dataspace);
  H5Dclose(dataset);


  /* write the attributes -- minimum and maximum values of the data */
  if (parallelIO ||
      (!parallelIO && (*myPE == MASTER_PE) && (*global_offset == 0))) {
    /* The logic ensures that creating and writing attributes work
       in both parallel and serial FLASH I/O modes.  In parallel all MPI
       processes create and write the attribute.  In serial the master
       MPI process creates and writes the attribute on the first
       entry to io_h5write_unknowns only.  Note: The io_h5_attribute_XXX
       functions open the dataset "record_label_new" before creating /
       writing attributes, so "record_label_new" should be closed before
       the call. */

    io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			   diskSize, record_label_new, minLabel);
    io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			  record_label_new, minLabel, varMin);

    io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			   diskSize, record_label_new, maxLabel);
    io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			  record_label_new, maxLabel, varMax);
  }
}
