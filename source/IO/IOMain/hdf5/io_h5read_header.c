/* This file contains the functions that read the data from the HDF5 file
 * The functions accept the PARAMESH data through arguments, since C cannot
 * handle common blocks 
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

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

/****if* source/IO/IOMain/hdf5/io_h5read_header_
 *
 * NAME
 *
 *  io_h5read_header_
 *
 *
 * SYNOPSIS
 *
 *  void io_h5read_header_(MyPE, file_identifier, unk_labels, outputSplitNum)
 *
 *  void io_h5read_header_(int *, hid_t *, int *, double *, double *, 
 *                            double *, int *)
 *
 *
 * DESCRIPTION
 *  Given an HDF5 file identifier (*file_identifier), read the FLASH header
 *  information from the file.  In this case the only array we really need is 
 *  the unknown variable names.
 *
 *
 * ARGUMENTS
 *
 *  MyPE              the local processor number
 *
 *  file_identifer    the HDF5 file handle (as returned from 
 *                    io_h5open_file_for_read_)
 *
 *  unk_labels        4 letter (5 if you count the null char) string given to 
 *                    unknown variable names
 *
 *  outputSplitNum    not implemented fully right now.  Needs more development.
 *                    Intended to allow users to split checkpoint file into x parts.
 *                    
 ****/

void FTOC(io_h5read_header)(int* MyPE,
                      hid_t* file_identifier,    /* file handle */
                      char unk_labels[][5],
                      int* outputSplitNum)      /*unknown labels */


{
  hid_t dataset;
  herr_t status;

  int file_version;

  hid_t sp_type;




}


