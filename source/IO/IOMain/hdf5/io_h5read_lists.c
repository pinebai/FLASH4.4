/* This file contains the functions that read the data from the HDF5 file
 * The functions accept the PARAMESH data through arguments, since C cannot
 * handle common blocks 
 */

 /* WARNING:  These functions side-effect their arguments.  The inputs *WILL* 
  * be overwritten with values appropriate to the checkpoint file being read
  * in as opposed to the value from the parfile/Config files.  Do not use the
  * values returned by these functions for writing.
  */

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

void io_h5readLists
  (hid_t* file_identifier,    /* file handle */
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
   char* dataset_names[])
{

  hsize_t dimens_1d, maxdimens_1d;

  hid_t dataspace, memspace, dataset;

  hid_t string_type;

  hid_t real_list_type;
  hid_t int_list_type;
  hid_t str_list_type;
  hid_t log_list_type;
  
  real_list_t *real_list;
  int_list_t *int_list;
  str_list_t *str_list;
  log_list_t *log_list;

  herr_t (*old_func)(void*) = NULL;
  void *old_client_data = NULL;
  herr_t status;

  int i;
  int errbufsize = 1024;
  char s[1024];

  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, MAX_STRING_LENGTH);


  H5Eget_auto(&old_func, &old_client_data); /* save old error handling here - KW */

  /*************************************
   * real runtime parameters / scalars *
   *************************************/

  /* temporarily turn off error printing */
  H5Eset_auto(NULL, NULL);
  /* create the name of the dataset */
  dataset = H5Dopen(*file_identifier, dataset_names[2]);

  /* Restore error handling */
  H5Eset_auto(old_func, old_client_data);
  
  if (dataset < 0){
    snprintf(s, errbufsize,"couldn't open the dataset %s\n", dataset_names[2]);
    Driver_abortFlashC(s);
  }

  dataspace = H5Dget_space(dataset);

  /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
  if (dimens_1d > *num_real) {
    snprintf(s,errbufsize,"Error reading: more than %d %s in checkpoint file!\n", *num_real, dataset_names[2]);
    Driver_abortFlashC(s);
  }
  *num_real = dimens_1d;

  /* malloc a pointer to a list of real_list_t's */
  real_list = (real_list_t *) malloc(dimens_1d * sizeof(real_list_t)); 

  /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
  memspace = H5Screate_simple(1, &dimens_1d, NULL);

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
  
  /* read the data into 'real_list' */
  status = H5Dread(dataset, real_list_type, memspace, dataspace,
                H5P_DEFAULT, real_list);


  if (status < 0) {
    snprintf(s,errbufsize,"Error reading %s from data file\n", dataset_names[2]);
    Driver_abortFlashC(s);
  }

  for (i = 0; i < dimens_1d; i++) {
    strncpy(real_names[i], real_list[i].name, MAX_STRING_LENGTH);
    real_values[i] = real_list[i].value;
  }

  free(real_list);
  H5Tclose(real_list_type);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);



  /****************************************
   * integer runtime parameters / scalars *
   ****************************************/
  
  /* temporarily turn off error printing */
  H5Eset_auto(NULL, NULL);
  /* create the name of the dataset */
  dataset = H5Dopen(*file_identifier, dataset_names[0]); 

  /* Restore error handling */
  H5Eset_auto(old_func, old_client_data);
    
  if (dataset < 0){
    snprintf(s,errbufsize,"couldn't open the dataset %s\n", dataset_names[0]);
    Driver_abortFlashC(s);
  }

  dataspace = H5Dget_space(dataset);

  /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
  if (dimens_1d > *num_int) {
    snprintf(s,errbufsize,"Error reading: more than %d %s in checkpoint file!\n", *num_int, dataset_names[0]);
    Driver_abortFlashC(s);
  }
  *num_int = dimens_1d;

  /* malloc a pointer to a list of int_list_t's */
  int_list = (int_list_t *) malloc(dimens_1d * sizeof(int_list_t)); 

  /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
  memspace = H5Screate_simple(1, &dimens_1d, NULL);

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
  
  /* read the data into 'int_list' */
  status = H5Dread(dataset, int_list_type, memspace, dataspace,
                H5P_DEFAULT, int_list);

  if (status < 0) {
    snprintf(s,errbufsize,"Error reading %s from data file\n", dataset_names[0]);
    Driver_abortFlashC(s);
  }
  
  for (i = 0; i < dimens_1d; i++) {
    strncpy(int_names[i], int_list[i].name, MAX_STRING_LENGTH);
    int_values[i] = int_list[i].value;
  }

  free(int_list);
  H5Tclose(int_list_type);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);



  /***************************************
   * string runtime parameters / scalars *
   ***************************************/

  /* temporarily turn off error printing */
  H5Eset_auto(NULL, NULL);
  /* create the name of the dataset */
  dataset = H5Dopen(*file_identifier, dataset_names[1]); 

  /* Restore error handling */
  H5Eset_auto(old_func, old_client_data);
  
  if (dataset < 0){
    snprintf(s,errbufsize,"couldn't open the dataset %s\n", dataset_names[1]);
    Driver_abortFlashC(s);
  }

  dataspace = H5Dget_space(dataset);

  /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
  if (dimens_1d > *num_str) {
    snprintf(s,errbufsize,"Error reading: more than %d %s in checkpoint file!\n", *num_str, dataset_names[1]);
    Driver_abortFlashC(s);
  }
  *num_str = dimens_1d;

  /* malloc a pointer to a list of str_list_t's */
  str_list = (str_list_t *) malloc(dimens_1d * sizeof(str_list_t)); 

  /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
  memspace = H5Screate_simple(1, &dimens_1d, NULL);

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
  
  /* read the data into 'str_list' */
  status = H5Dread(dataset, str_list_type, memspace, dataspace,
                H5P_DEFAULT, str_list);


  if (status < 0) {
    snprintf(s,errbufsize,"Error reading %s from data file\n", dataset_names[1]);
    Driver_abortFlashC(s);
  }
  
  for (i = 0; i < dimens_1d; i++) {
    strncpy(str_names[i], str_list[i].name, MAX_STRING_LENGTH);
    strncpy(str_values[i], str_list[i].value, MAX_STRING_LENGTH);
  }

  free(str_list);
  H5Tclose(str_list_type);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);



  /****************************************
   * logical runtime parameters / scalars *
   ****************************************/

  /* temporarily turn off error printing */
  H5Eset_auto(NULL, NULL);
  /* create the name of the dataset */
  dataset = H5Dopen(*file_identifier, dataset_names[3]); 

  /* Restore error handling */
  H5Eset_auto(old_func, old_client_data);
  
  if (dataset < 0){
    snprintf(s,errbufsize,"couldn't open the dataset %s\n", dataset_names[3]);
    Driver_abortFlashC(s);
  }

  dataspace = H5Dget_space(dataset);

  /* read extent of 'dataspace' (i.e. # of name/value pairs) into 'dimens_1d' */
  H5Sget_simple_extent_dims(dataspace, &dimens_1d, &maxdimens_1d);
  if (dimens_1d > *num_log) {
    snprintf(s,errbufsize,"Error reading: more than %d %s in checkpoint file!\n", *num_log, dataset_names[3]);
    Driver_abortFlashC(s);
  }
  *num_log = dimens_1d;

  /* malloc a pointer to a list of log_list_t's */
  log_list = (log_list_t *) malloc(dimens_1d * sizeof(log_list_t)); 

  /* create a new simple dataspace of 1 dimension and size of 'dimens_1d' */
  memspace = H5Screate_simple(1, &dimens_1d, NULL);

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
  
  /* read the data into 'log_list' */
  status = H5Dread(dataset, log_list_type, memspace, dataspace,
                H5P_DEFAULT, log_list);


  if (status < 0) {
    snprintf(s,errbufsize,"Error reading %s from data file\n", dataset_names[3]);
    Driver_abortFlashC(s);
  }
  
  for (i = 0; i < dimens_1d; i++) {
    strncpy(log_names[i], log_list[i].name, MAX_STRING_LENGTH);
    log_values[i] = log_list[i].value;
  }

  free(log_list);
  H5Tclose(log_list_type);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  H5Tclose(string_type);
}





void FTOC(io_h5read_runtime_parameters)
     (hid_t* file_identifier,    /* file handle */
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
      int log_values[])
{

  char* dataset_names[] = {"integer runtime parameters",
                           "string runtime parameters",
                           "real runtime parameters",
                           "logical runtime parameters"};


  io_h5readLists
     (file_identifier,
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
      dataset_names);
}




void FTOC(io_h5read_scalars)
     (hid_t* file_identifier,    /* file handle */
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
      int log_values[])
{

  char* dataset_names[] = {"integer scalars",
                           "string scalars",
                           "real scalars",
                           "logical scalars"};


  io_h5readLists
     (file_identifier,
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
      dataset_names);
}
