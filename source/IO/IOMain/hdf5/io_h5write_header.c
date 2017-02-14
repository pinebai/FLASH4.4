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
#include <assert.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/io_h5write_header_
 *
 * NAME
 *
 *  io_h5write_header_
 *
 *
 * SYNOPSIS
 *
 *
 *  void FTOC(io_h5write_header)(int* MyPE,
 *                     int* nvar_out,
 *                     hid_t* file_identifier,
 *                     char geometry_name[],    
 *                     char unk_labels[][4],  
 *                     char setup_call[],  
 *                     char file_creation_time[],
 *                     char flash_version[],     
 *                     char build_date[],  
 *                     char build_dir[],   
 *                     char build_machine[],
 *                     char cflags[],
 *                     char fflags[],
 *                     char setup_time_stamp[],
 *                     char build_time_stamp[],
 *                     int *pFileFormat,
 *                     int* io_outputSplitNum)
 *
 *
 *
 * DESCRIPTION
 *
 *  Write the header information describing the FLASH dataset into the HDF5 
 *  file.
 *
 *
 * ARGUMENTS
 *
 *  MyPE                  processor number
 *
 *  nvar_out              the number of solution variables to store
 *
 *  file_identifier       the file handle as return by the HDF5 library
 *
 *  unk_labels            the 4-character descriptive name for each of the 
 *                        solution variables we are storing.
 *
 *  setup_call            a string giving the syntax of the setup command used
 *                        when the current executable was built
 *
 *
 *  file_creation_time    a string giving the current time/date
 *
 *  flash_version         a string giving the version of the FLASH source used
 *
 *  build_date            a string giving the date that the executable 
 *                        was built
 *
 *  build_dir             a string giving the path of the build directory
 *
 *  build_machine         a string giving the uname information of the build 
 *                        machine
 *
 *  fflags                 f compiler flags
 *
 *  cflags                 c compiler flags
 *
 *  setup_time_stamp       time when simulation was setup
 *
 *  build_time_stamp       time when simulation was built
 *
 *  pFileFormat            Flash file format
 *
 *  io_outputSplitNum      not fully implemented yet.  Intended to allow user to
 *                         split checkpoint files into x parts
 *
 ****/

void FTOC(io_h5write_header)(int* MyPE,
                       int* nvar_out,             /* num vars to store */
                       hid_t* file_identifier,    /* file handle */
                       char geometry_name[],    
                       char unk_labels[][4],      /* unknown labels */
                       char setup_call[],   /* syntax of setup */
                       char file_creation_time[], /* time and date stamp */
                       char flash_version[],      /* FLASH version num */
                       char build_date[],   /* date of compilation */
                       char build_dir[],    /* path of compilation */
                       char build_machine[], /* machine compiled on */
                       char cflags[],
                       char fflags[],
                       char setup_time_stamp[],
                       char build_time_stamp[],
		       int *pFileFormat,
                       int* io_outputSplitNum)

{


  hid_t dataspace, dataset;

  herr_t status;

  int rank;
  hsize_t dimens_1d, dimens_2d[2];

  int string_size;

  hid_t string_type, sp_type, string_type_setup, string_type_max_str_len, si_type;


  sim_info_t sim_info;

  int file_version;

  hid_t attribute, attribute_space;

  hid_t group_id;
  int aerror, i,j;
  


  /*-----------------------------------------------------------------------
    create a compound datatype to hold the simulation info
    These are things like file creation time, file formate version
    flash version, build date
    build dir, build machine, setup call */
  
  


  string_type_max_str_len = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_max_str_len, MAX_STRING_LENGTH);



  string_type_setup = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type_setup, 400);

  sim_info.file_format_version = *pFileFormat;
#ifndef FLASH_IO_EXPERIMENTAL
  assert(*pFileFormat == FILE_FORMAT_VERSION);
#endif


  strncpy(sim_info.setup_call, setup_call, 400);

  strncpy(sim_info.file_creation_time, file_creation_time, MAX_STRING_LENGTH);

  strncpy(sim_info.flash_version, flash_version, MAX_STRING_LENGTH);

  strncpy(sim_info.build_date, build_date, MAX_STRING_LENGTH);

  strncpy(sim_info.build_dir, build_dir, MAX_STRING_LENGTH);

  strncpy(sim_info.build_machine, build_machine, MAX_STRING_LENGTH);

  strncpy(sim_info.cflags, cflags, 400);

  strncpy(sim_info.fflags, fflags, 400);

  strncpy(sim_info.setup_time_stamp, setup_time_stamp, MAX_STRING_LENGTH);

  strncpy(sim_info.build_time_stamp, build_time_stamp, MAX_STRING_LENGTH);



  
  rank = 1;
  dimens_1d = 1;

  
  dataspace = H5Screate_simple(rank, &dimens_1d, NULL);

  si_type = H5Tcreate(H5T_COMPOUND, sizeof(sim_info_t));

  
#undef offsetof
#define offsetof(st, m) ((size_t) ( (char *)&((st *)(0))->m - (char *)0 ))

  H5Tinsert(si_type, 
          "file format version", 
          offsetof(sim_info_t, file_format_version),
          H5T_NATIVE_INT);
  
  H5Tinsert(si_type, 
          "setup call", 
          HOFFSET(sim_info_t, setup_call),
          string_type_setup);
  
  H5Tinsert(si_type, 
          "file creation time", 
          HOFFSET(sim_info_t, file_creation_time),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "flash version", 
          HOFFSET(sim_info_t, flash_version),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "build date", 
          HOFFSET(sim_info_t, build_date),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "build dir", 
          HOFFSET(sim_info_t, build_dir),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "build machine", 
          HOFFSET(sim_info_t, build_machine),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "cflags", 
          HOFFSET(sim_info_t, cflags),
          string_type_setup);
  

  H5Tinsert(si_type, 
          "fflags", 
          HOFFSET(sim_info_t, fflags),
          string_type_setup);
  

  H5Tinsert(si_type, 
          "setup time stamp", 
          HOFFSET(sim_info_t,setup_time_stamp),
          string_type_max_str_len);
  
  H5Tinsert(si_type, 
          "build time stamp", 
          HOFFSET(sim_info_t, build_time_stamp),
          string_type_max_str_len);
  
  
  
  dataset = H5Dcreate(*file_identifier, "sim info", si_type,
                  dataspace, H5P_DEFAULT);

  

  if ((*MyPE == 0) || (*io_outputSplitNum != 1)) {
    status = H5Dwrite(dataset, si_type, H5S_ALL, H5S_ALL,
                  H5P_DEFAULT, &sim_info);
  }
  

  H5Tclose(string_type_setup);
  H5Tclose(string_type_max_str_len);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(si_type);
  
  

   /*------------------------------------------------------------------------- */


#ifndef FLASH_IO_EXPERIMENTAL
  /* unknown names */

  if(*nvar_out > 0){

    /* first we have to lowercase all of the names so that these are 
       consistent with the unkNames.*/
    for (i = 0; i < *nvar_out; ++i){
      for (j=0; j<4; ++j){
	if( unk_labels[i][j] >= 'A' && unk_labels[i][j] <= 'Z')
	  unk_labels[i][j] += 32;
      }
    }

    rank = 2;
    dimens_2d[0] = (hsize_t) *nvar_out;
    dimens_2d[1] = 1;
    
    
    /* manually set the string size */
    string_size = 4;
    
    /* setup the datatype for this string length */
    string_type = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(string_type, string_size);
    
#ifdef DEBUG_IO
    printf("MyPE = %d, string type = %d, set size status = %d\n",
         *MyPE, (int) string_type, (int) status);
#endif
    
    dataspace = H5Screate_simple(rank, dimens_2d, NULL);
    dataset   = H5Dcreate(*file_identifier, "unknown names", 
                    string_type, dataspace, H5P_DEFAULT);
    
    if ((*MyPE == 0) || (*io_outputSplitNum != 1)) {
      status    = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, 
                     H5P_DEFAULT, unk_labels);
    }
    H5Tclose(string_type);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
#endif

}


