/* This file contains the functions that write the data to the HDF5 file
 * The functions accept the PARAMESH data through arguments, since C cannot
 * handle common blocks 
 */

/* if the following flag is defined, status checking on the writes will
   be performed, and the results will be output to stdout */

/* #define DEBUG_IO */

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
/* This function writes lists of the runtime parameters and the scalars.
   Both are stored in the same way under the hood. The runtime parameters
   and scalars are each put into separate lists based on type (int, str,
   real or logical) These values are then written out separately and stored
   in the checkpoint file as "int runtime parameters / scalars", "real 
   runtime parameters / scalars", etc. */

void io_h5writeLists
  (int* myPE, 
   hid_t* file_identifier,    /* file handle */
   int* num_real, /* number of real values */
   char real_names[][MAX_STRING_LENGTH],
   double real_values[], 
   int* num_int, 
   char int_names[][MAX_STRING_LENGTH], 
   int int_values[],
   int* num_str,
   char str_names[][MAX_STRING_LENGTH],
   char str_values[][MAX_STRING_LENGTH],
   int* num_log, 
   char log_names[][MAX_STRING_LENGTH],
   int log_values[],
   char* dataset_names[],
   int* io_outputSplitNum)
     
{
  hid_t dataspace, dataset;
  herr_t status;

  hsize_t dimens_1d;  

  hid_t string_type;
  
  hid_t real_list_type;
  hid_t int_list_type;
  hid_t str_list_type;
  hid_t log_list_type;
  
  real_list_t *real_list;
  int_list_t *int_list;
  str_list_t *str_list;
  log_list_t *log_list;
  
  int error_count = 0;
  int i;

  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, MAX_STRING_LENGTH);


  /*************************************
   * real runtime parameters / scalars *
   *************************************/

  if (*num_real > 0) {

    /* malloc a pointer to a list of real_list_t's */
    real_list = (real_list_t *) malloc(*num_real * sizeof(real_list_t)); 
  
    /* copy names and values into that list */  
    for (i = 0; i < *num_real; i++) {
      strncpy(real_list[i].name, real_names[i], MAX_STRING_LENGTH);
      real_list[i].value = real_values[i];
    }

    dimens_1d = *num_real;

    /* create a new simple dataspace of 1 dimen. and size of 'num_real' */    
    dataspace = H5Screate_simple(1, &dimens_1d, NULL);

    /* create an empty vessel sized to hold one real_list_t's worth of data */    
    real_list_type = H5Tcreate(H5T_COMPOUND, sizeof(real_list_t));

#undef offsetof
#define offsetof(st, m) ((size_t) ( (char *)&((st *)(0))->m - (char *)0 ))

    /* subdivide the empty vessel into its component sections (name and value) */    
    H5Tinsert(real_list_type, 
            "name", 
            HOFFSET(real_list_t, name),
            string_type);
    
    H5Tinsert(real_list_type, 
            "value", 
            HOFFSET(real_list_t, value),
            H5T_NATIVE_DOUBLE);

    dataset = H5Dcreate(*file_identifier, dataset_names[2], real_list_type,
                  dataspace, H5P_DEFAULT);

    if ((*myPE == 0) || (*io_outputSplitNum != 1)) {
      /* write the data into the checkpoint file */    
      status = H5Dwrite(dataset, real_list_type, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, real_list);
      
      if (status < 0) {
      printf("Error writing %s to checkpoint file\n", dataset_names[2]);
      error_count++;
      }
    }    

    free(real_list);
    H5Tclose(real_list_type);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }



  /****************************************
   * integer runtime parameters / scalars *
   ****************************************/
  
  if (*num_int > 0) {

    /* malloc a pointer to a list of int_list_t's */
    int_list = (int_list_t *) malloc(*num_int * sizeof(int_list_t)); 
  
    /* copy names and values into that list */  
    for (i = 0; i < *num_int; i++) {
      strncpy(int_list[i].name, int_names[i], MAX_STRING_LENGTH);
      int_list[i].value = int_values[i];
    }
 
    dimens_1d = *num_int;

    /* create a new simple dataspace of 1 dimen. and size of 'num_int' */    
    dataspace = H5Screate_simple(1, &dimens_1d, NULL);

    /* create an empty vessel sized to hold one int_list_t's worth of data */    
    int_list_type = H5Tcreate(H5T_COMPOUND, sizeof(int_list_t));

    /* subdivide the empty vessel into its component sections (name and value) */    
    H5Tinsert(int_list_type, 
            "name", 
            HOFFSET(int_list_t, name),
            string_type);
    
    H5Tinsert(int_list_type, 
            "value", 
            HOFFSET(int_list_t, value),
            H5T_NATIVE_INT);

    dataset = H5Dcreate(*file_identifier, dataset_names[0], int_list_type,
                  dataspace, H5P_DEFAULT);

    if ((*myPE == 0) || (*io_outputSplitNum != 1)) {
      /* write the data into the checkpoint file */    
      status = H5Dwrite(dataset, int_list_type, H5S_ALL, H5S_ALL,
                        H5P_DEFAULT, int_list);
      
      if (status < 0) {
        printf("Error writing %s to checkpoint file\n", dataset_names[0]);
        error_count++;
      }
    }
    
    free(int_list);
    H5Tclose(int_list_type);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }



  /***************************************
   * string runtime parameters / scalars *
   ***************************************/

  if (*num_str > 0) {

    /* malloc a pointer to a list of str_list_t's */
    str_list = (str_list_t *) malloc(*num_str * sizeof(str_list_t)); 
  
    /* copy names and values into that list */  
    for (i = 0; i < *num_str; i++) {
      strncpy(str_list[i].name, str_names[i], MAX_STRING_LENGTH);
      strncpy(str_list[i].value, str_values[i], MAX_STRING_LENGTH);
    }

    dimens_1d = *num_str;

    /* create a new simple dataspace of 1 dimen. and size of 'num_str' */    
    dataspace = H5Screate_simple(1, &dimens_1d, NULL);

    /* create an empty vessel sized to hold one str_list_t's worth of data */    
    str_list_type = H5Tcreate(H5T_COMPOUND, sizeof(str_list_t));

    /* subdivide the empty vessel into its component sections (name and value) */    
    H5Tinsert(str_list_type, 
            "name", 
            HOFFSET(str_list_t, name),
            string_type);
    
    H5Tinsert(str_list_type, 
            "value", 
            HOFFSET(str_list_t, value),
            string_type);

    dataset = H5Dcreate(*file_identifier, dataset_names[1], str_list_type,
                  dataspace, H5P_DEFAULT);

    if ((*myPE == 0) || (*io_outputSplitNum != 1))  {
      /* write the data into the checkpoint file */    
      status = H5Dwrite(dataset, str_list_type, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, str_list);
      
      if (status < 0) {
      printf("Error writing %s to checkpoint file\n", dataset_names[1]);
      error_count++;
      }
    }

    free(str_list);
    H5Tclose(str_list_type);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }



  /****************************************
   * logical runtime parameters / scalars *
   ****************************************/

  if (*num_log > 0) {

    /* malloc a pointer to a list of log_list_t's */
    log_list = (log_list_t *) malloc(*num_log * sizeof(log_list_t)); 
  
    /* copy names and values into that list */  
    for (i = 0; i < *num_log; i++) {
      strncpy(log_list[i].name, log_names[i], MAX_STRING_LENGTH);
      log_list[i].value = log_values[i];
    }

    dimens_1d = *num_log;

    /* create a new simple dataspace of 1 dimen. and size of 'num_log' */    
    dataspace = H5Screate_simple(1, &dimens_1d, NULL);

    /* create an empty vessel sized to hold one log_list_t's worth of data */    
    log_list_type = H5Tcreate(H5T_COMPOUND, sizeof(log_list_t));

    /* subdivide the empty vessel into its component sections (name and value) */    
    H5Tinsert(log_list_type, 
            "name", 
            HOFFSET(log_list_t, name),
            string_type);
    
    H5Tinsert(log_list_type, 
            "value", 
            HOFFSET(log_list_t, value),
            H5T_NATIVE_INT);

    dataset = H5Dcreate(*file_identifier, dataset_names[3], log_list_type,
                  dataspace, H5P_DEFAULT);

    if ((*myPE == 0) || (*io_outputSplitNum != 1)) {
      /* write the data into the checkpoint file */    
      status = H5Dwrite(dataset, log_list_type, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, log_list);
      
      if (status < 0) {
      printf("Error writing %s to checkpoint file\n", dataset_names[3]);
      error_count++;
      }
    }
    free(log_list);
    H5Tclose(log_list_type);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }
  H5Tclose(string_type);
}


void FTOC(io_h5write_runtime_parameters)
     (int* myPE, 
      hid_t* file_identifier,    /* file handle */
      int* num_real,
      char real_names[][MAX_STRING_LENGTH],
      double real_values[], 
      int* num_int, 
      char int_names[][MAX_STRING_LENGTH], 
      int int_values[],
      int* num_str,
      char str_names[][MAX_STRING_LENGTH],
      char str_values[][MAX_STRING_LENGTH],
      int* num_log, 
      char log_names[][MAX_STRING_LENGTH],
      int log_values[],
      int* io_outputSplitNum)

{

  char* dataset_names[] = {"integer runtime parameters",
                           "string runtime parameters",
                           "real runtime parameters",
                           "logical runtime parameters"};


  io_h5writeLists
     (myPE, 
      file_identifier,
      num_real,
      real_names,
      real_values, 
      num_int, 
      int_names, 
      int_values,
      num_str,
      str_names,
      str_values,
      num_log, 
      log_names,
      log_values,
      dataset_names,
      io_outputSplitNum);
}


void FTOC(io_h5write_scalars)
     (int* myPE, 
      hid_t* file_identifier,    /* file handle */
      int* num_real,
      char real_names[][MAX_STRING_LENGTH],
      double real_values[], 
      int* num_int, 
      char int_names[][MAX_STRING_LENGTH], 
      int int_values[],
      int* num_str,
      char str_names[][MAX_STRING_LENGTH],
      char str_values[][MAX_STRING_LENGTH],
      int* num_log, 
      char log_names[][MAX_STRING_LENGTH],
      int log_values[],
      int* io_outputSplitNum)
{

  char* dataset_names[] = {"integer scalars",
                           "string scalars",
                           "real scalars",
                           "logical scalars"};


  io_h5writeLists
     (myPE, 
      file_identifier,
      num_real,
      real_names,
      real_values, 
      num_int, 
      int_names, 
      int_values,
      num_str,
      str_names,
      str_values,
      num_log, 
      log_names,
      log_values,
      dataset_names,
      io_outputSplitNum);
}
