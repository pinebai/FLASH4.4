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
#include <pnetcdf.h>
#include "mangle_names.h"
#include <mpi.h>
#include "Flash.h"
#include "constants.h"

#define FILE_FORMAT_VERSION 9

int Driver_abortFlashC(char* message);



void FTOC(io_ncmpi_write_header)(int* MyPE,
			 int* pFile_version,
                         int* nvar_out, 
                         int nzones_block[3],  
                         char labels[][4],
                         int* file_identifier,
                         int* total_blocks,
                         int* ngid,
                         int* numRealParms,
                         char realParmNames[][MAX_STRING_LENGTH], 
                         double realParmValues[],
                         int* numIntParms,
                         char intParmNames[][MAX_STRING_LENGTH], 
                         int intParmValues[],
                         int* numStrParms,
                         char strParmNames[][MAX_STRING_LENGTH], 
                         char strParmValues[][MAX_STRING_LENGTH], 
                         int* numLogParms,
                         char logParmNames[][MAX_STRING_LENGTH], 
                         int logToIntParmValues[],
                         int* numRealScalars,
                         char realScalarNames[][MAX_STRING_LENGTH], 
                         double realScalarValues[],
                         int* numIntScalars,
                         char intScalarNames[][MAX_STRING_LENGTH], 
                         int intScalarValues[],
                         int* numStrScalars,
                         char strScalarNames[][MAX_STRING_LENGTH], 
                         char strScalarValues[][MAX_STRING_LENGTH], 
                         int* numLogScalars,
                         char logScalarNames[][MAX_STRING_LENGTH], 
                         int logToIntScalarValues[],
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
			 int* doublePrecisionInt,
			 int* metadataDP,
			 int* gridQuantityDP,
			 int* globalNumParticles) 


{

  int i,k;
  int ncid, status = 0;
  int rank;
  int dim_tot_blocks, dim_nxb, dim_nyb, dim_nzb, dim_NDIM, dim_NGID, dim_2, dim_1;
  int dim_MAX_STR_LEN, dim_MDIM;
  int dimids[5]; 
  int two = 2;
  int one = 1;
  char record_label[5];
  int varid[500];
  int id = 0;
  /*char labelstring[25];*/
  char str[4]; 
  char mkeys[12];
  char *p;
  int long_len = 400;
  int file_version;
  double varMinDummy, varMaxDummy;
  float varMinDummyFloat, varMaxDummyFloat;


  

  ncid = *file_identifier;
  

  
  status = ncmpi_def_dim(ncid, "dim_tot_blocks", (MPI_Offset)(*total_blocks), &dim_tot_blocks);
  if (status < 0){
    printf("Error: Unable to define dim_tot_blocks\n");
    Driver_abortFlashC("Error: Unable to define dim_tot_blocks\n");
  }
 
   /*Added to properly handle face-centered variable data and potentially scratch data*/
   status = ncmpi_def_dim(ncid, "dim_nxb_long", (MPI_Offset)nzones_block[0]+1, &dim_nxb);
  if (status < 0){
    printf("Error: Unable to define dim_nxb_long\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nyb_long", (MPI_Offset)nzones_block[1]+1, &dim_nyb);
  if (status < 0){
    printf("Error: Unable to define dim_nyb\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nzb_long", (MPI_Offset)nzones_block[2]+1, &dim_nzb);
  if (status < 0){
    printf("Error: Unable to define dim_nzb\n");
  }


  status = ncmpi_def_dim(ncid, "dim_nxb", (MPI_Offset)nzones_block[0], &dim_nxb);
  if (status < 0){
    printf("Error: Unable to define dim_nxb\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nyb", (MPI_Offset)nzones_block[1], &dim_nyb);
  if (status < 0){
    printf("Error: Unable to define dim_nyb\n");
  }
  status = ncmpi_def_dim(ncid, "dim_nzb", (MPI_Offset)nzones_block[2], &dim_nzb);
  if (status < 0){
    printf("Error: Unable to define dim_nzb\n");
  }



  status = ncmpi_def_dim(ncid, "dim_NGID", (MPI_Offset)*ngid, &dim_NGID);
  if (status < 0){
    printf("Error: Unable to define dim_NGID\n");
    Driver_abortFlashC("Error: Unable to define dim_NGID\n");
  }

  status = ncmpi_def_dim(ncid, "dim_NDIM", (MPI_Offset)NDIM, &dim_NDIM);
  if (status < 0){
    printf("Error: Unable to define dim_NDIM\n");
    Driver_abortFlashC("Error: Unable to define dim_NDIM\n");
  }

  status = ncmpi_def_dim(ncid, "dim_2", (MPI_Offset)two, &dim_2);
  if (status < 0){
    printf("Error: Unable to define dim_2\n");
    Driver_abortFlashC("Error: Unable to define dim_2\n");
  }


  status = ncmpi_def_dim(ncid, "dim_1", (MPI_Offset)one, &dim_1);
  if (status < 0){
    printf("Error: Unable to define dim_1\n");
    Driver_abortFlashC("Error: Unable to define dim_1\n");
  }


  status = ncmpi_def_dim(ncid, "dim_MAX_STR_LEN", (MPI_Offset)MAX_STRING_LENGTH, &dim_MAX_STR_LEN);
  if (status < 0){
    printf("Error: Unable to define dim_MAX_STR_LEN\n");
    Driver_abortFlashC("Error: Unable to define dim_MAX_STR_LEN\n");
  }

  status = ncmpi_def_dim(ncid, "dim_MDIM", (MPI_Offset)MDIM, &dim_MDIM);
   if (status < 0){
    printf("Error: Unable to define dim_MDIM\n");
    Driver_abortFlashC("Error: Unable to define dim_MDIM\n");
  } 







  /*****************************************************************************/
  /* Define the vars */


  /* define holder var for runtime parameters 
     store all the runtime parameters together is some 'non-global'
     attribute.  Make a variable runtime_parameters and store the
     actually parameters as attributes of this variable */
  /*0 not used right now*/
  rank = 1;
  dimids[0] = dim_1;
  status = ncmpi_def_var (ncid, "runtime_parameters", NC_INT, rank, dimids, &varid[id]);
  /*printf("varid[%d] = %d\n", id, varid[id]);*/
  if (status < 0){
    printf("Error: Unable to define lrefine\n");
    Driver_abortFlashC("Error: Unable to define lrefine\n");
  }


  /* define holder var for scalars 
     same idea for scalar values.  Group them together and make
     the actual scalars attributes of this scalar variable */
  /* 1 */
  rank = 1;
  id ++;
  status = ncmpi_def_var (ncid, "scalars", NC_INT, rank, dimids, &varid[id]);
  /*printf("varid[%d] = %d\n", id, varid[id]);*/
  if (status < 0){
    printf("Error: Unable to define lrefine\n");
    Driver_abortFlashC("Error: Unable to define lrefine\n");
  }


  dimids[0] = dim_tot_blocks;

  /* 2 */
  /* define var for refinement level */
  rank = 1;
  id ++;
  status = ncmpi_def_var (ncid, "lrefine", NC_INT, rank, dimids, &varid[id]);
  /*printf("varid[%d] = %d\n", id, varid[id]);*/
  if (status < 0){
    printf("Error: Unable to define lrefine\n");
    Driver_abortFlashC("Error: Unable to define lrefine\n");
  }


#ifdef FLASH_GRID_PARAMESH3OR4
  /* define var for refinement level */
  rank = 1;
  id ++;
  status = ncmpi_def_var (ncid, "which_child", NC_INT, rank, dimids, &varid[id]);
  if (status < 0){
    printf("Error: Unable to define which_child\n");
    Driver_abortFlashC("Error: Unable to define which_child\n");
  }

  /* define var for gsurr_blks */
  rank = 5;
  id ++;
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_nzb;
  dimids[2] = dim_nyb;
  dimids[3] = dim_nxb;
  dimids[4] = dim_2;
  status = ncmpi_def_var (ncid, "gsurr_blks", NC_INT, rank, dimids, &varid[id]);
  if (status < 0){
    printf("Error: Unable to define gsurr_blks\n");
    Driver_abortFlashC("Error: Unable to define gsurr_blks\n");
  }
#endif



  /* define var for nodetype */
  rank = 1;
  id ++;
  status = ncmpi_def_var (ncid, "nodetype", NC_INT, rank, dimids, &varid[id]);
  if (status < 0){
    printf("Error: Unable to define nodetype\n");
    Driver_abortFlashC("Error: Unable to define nodetype\n");
  }



  /* define var for global id */
  rank = 2;
  id ++;
  dimids[1] = dim_NGID;
  status = ncmpi_def_var (ncid, "gid", NC_INT, rank, dimids, &varid[id]);
  if (status < 0){
    printf("Error: Unable to define gid\n");
    Driver_abortFlashC("Error: Unable to define gid\n");
  }


  /* define var for coordinates */
  rank = 2;
  id ++;
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_MDIM;
  /* if writing a plotfile values stored in single precision */
  if(*doublePrecisionInt || *metadataDP) {
    status = ncmpi_def_var (ncid, "coordinates", NC_DOUBLE, rank, dimids, &varid[id]);
  }else{
    status = ncmpi_def_var (ncid, "coordinates", NC_FLOAT, rank, dimids, &varid[id]);
  }

  if (status < 0){
    printf("Error: Unable to define coordinates\n");
    Driver_abortFlashC("Error: Unable to define coordinates\n");
  }  


  /* define var for blocksize */
  rank = 2;
  id ++;
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_MDIM;
  if(*doublePrecisionInt || *metadataDP) {
    status = ncmpi_def_var (ncid, "blocksize", NC_DOUBLE, rank, dimids, &varid[id]);
  }else{
    status = ncmpi_def_var (ncid, "blocksize", NC_FLOAT, rank, dimids, &varid[id]);
  }
  if (status < 0){
    printf("Error: Unable to define blocksize\n");
    Driver_abortFlashC("Error: Unable to define blocksize\n");
  }



  /* define var for bndbox */
  rank = 3;
  id ++;
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_MDIM;
  dimids[2] = dim_2;
  if(*doublePrecisionInt || *metadataDP) {
    status = ncmpi_def_var (ncid, "bndbox", NC_DOUBLE, rank, dimids, &varid[id]);
  }else{
    status = ncmpi_def_var (ncid, "bndbox", NC_FLOAT, rank, dimids, &varid[id]);
  }
  if (status < 0){
    printf("Error: Unable to define bndbox\n");
    Driver_abortFlashC("Error: Unable to define bndbox\n");
  }


  /* define var for processor number */
  rank = 1;
  id ++;
  dimids[0] = dim_tot_blocks;
  status = ncmpi_def_var (ncid, "processor_number", NC_INT, rank, dimids, &varid[id]);
  if (status < 0){
    printf("Error: Unable to define proc_number\n");
    Driver_abortFlashC("Error: Unable to define proc_number\n");
  }


 
#ifndef FLASH_IO_EXPERIMENTAL
  /* define var for unknown array */
  
  rank = 4;
  dimids[0] = dim_tot_blocks;
  dimids[1] = dim_nzb;
  dimids[2] = dim_nyb;
  dimids[3] = dim_nxb;
  

  varMinDummy = 999.9;
  varMaxDummy = 999.9;
  varMinDummyFloat = 999.;
  varMaxDummyFloat = 999.;


  
  for (i=0; i<*nvar_out; i++){
    strncpy(record_label, labels[i], 4);
    for(k=0; k<4; k++) {
      if(record_label[k] == ' ') {
      record_label[k] = '_';
      }
    }
    *(record_label+4) = '\0';

    id ++;
    if((*doublePrecisionInt == 1) || *gridQuantityDP) {
      status = ncmpi_def_var (ncid, record_label, NC_DOUBLE, rank, dimids, &varid[id]);

      if (status < 0){
	printf("Error: Unable to define %s\n", record_label);
	Driver_abortFlashC("Error: Unable to define var record_label\n");
      }
      

      /* defining min and max attributes for each  unknown here in the header */
      /* filling the values with dummy varMin and varMax */
      /* we actually write the correct values in io_ncmpi_write_header by reentering define mode */
      /* it was recommended to do it this way by the pnetcdf guys when we saw a performance */
      /* hit reentering define mode. */
	 
      status = ncmpi_put_att_double(ncid, varid[id], "minimum", NC_DOUBLE, 1, &varMinDummy);
      if (status < 0){
	printf("Error: ncmpi_put_att_double minimum");
	Driver_abortFlashC("Error: ncmpi_put_att_double, minimum\n");
      }

      status = ncmpi_put_att_double(ncid, varid[id], "maximum", NC_DOUBLE, 1, &varMaxDummy);
      if (status < 0){
	printf("Error: ncmpi_put_att_double maximum");
	Driver_abortFlashC("Error: ncmpi_put_att_double, maximum\n");
      }



    }else{
      status = ncmpi_def_var (ncid, record_label, NC_FLOAT, rank, dimids, &varid[id]);
       
    if (status < 0){
      printf("Error: Unable to define %s\n", record_label);
      Driver_abortFlashC("Error: Unable to define var record_label\n");
    }

    /* See note above about putting attributes for unknowns in first define mode */
    /* Correct value of varMinDummyFloat and varMaxDummyFloat written in */
    /* io_ncmpi_write_unknowns_sp */
    status = ncmpi_put_att_float(ncid, varid[id], "minimum", NC_FLOAT, 1, &varMinDummyFloat);
    if (status < 0){
      printf("Error: io_ncmpi_write_header ncmpi_put_att_double minimum");
      Driver_abortFlashC("Error: io_ncmpi_write_header ncmpi_put_att_double, minimum\n");
    }
    
    status = ncmpi_put_att_float(ncid, varid[id], "maximum", NC_FLOAT, 1, &varMaxDummyFloat);
    if (status < 0){
      printf("Error: io_ncmpi_write_header ncmpi_put_att_double maximum");
      Driver_abortFlashC("Error:io_ncmpi_write_header ncmpi_put_att_double, maximum\n");
    }
    }

  }
  



  /* Write out modules used 
     ncmpi_put_att_int(ncid, NC_GLOBAL, "number_modules", NC_INT, 1, num_modules);
     for(i=0; i<*num_modules; i++) {
     sprintf(str, "%d", i);
     strcpy(mkeys, "module_");
     strcat(mkeys,str);
     p = strtok(module_names[i], "    ");
     string_size=strlen(p);
     status = ncmpi_put_att_text(ncid, NC_GLOBAL, mkeys, string_size, p);
     if (status < 0){
     printf("Error: Unable to write string_size\n");
     Driver_abortFlashC("Error: Unable to write string_size\n");
     }
     }
  */
  
  /*Write out unknown names*/ 
  
  ncmpi_put_att_int(ncid, NC_GLOBAL, "unknown_names", NC_INT, 1, nvar_out);
  for(i=0; i<*nvar_out; i++) {
    sprintf(str, "%d", i);
    strcpy(mkeys, "unknown_");
    strcat(mkeys,str);
    strncpy(record_label, labels[i], 4);
    for(k=0; k<=4; k++) {
      if(record_label[k] == ' ') {
      record_label[k] = '_';
      }
    }
    *(record_label+4) = '\0';
    status = ncmpi_put_att_text(ncid, NC_GLOBAL, mkeys, 4, record_label);
    if (status < 0){
      printf("Error: ncmpi_put_att_text\n");
      Driver_abortFlashC("Error: ncmpi_put_att_text\n");
    }
  }


  /*
    for(i=0; i<*num_modules; i++) {
    sprintf(str, "%d", i);
    p = strtok(module_names[i], "    ");
    string_size=strlen(p);
    status = ncmpi_put_att_text(ncid, NC_GLOBAL, str, string_size, p);
    if (status < 0){
    printf("Error: ncmpi_put_att_text\n");
    Driver_abortFlash("Error: ncmpi_put_att_text\n");
    }
    }

    write unknown extremas as attributes
    maximums
    for(i=0; i <*nvar_out; i++) {
    strncpy(record_label, labels[i], 4);
            for(k=0; k<=4; k++) {
                      if(record_label[k] == ' ') {
                                   record_label[k] = '_';
                                 }
                  }
            *(record_label+4) = '\0';
            p = strcat(record_label,max);
            status = ncmpi_put_att_double(ncid, NC_GLOBAL, p, NC_DOUBLE, 1, &varMax[i]);
            if (status < 0){
              printf("Error: ncmpi_put_att_text, varMax\n");
              Driver_abortFlash("Error: ncmpi_put_att_text, varMax\n");
            }
            }
      minimums
      for(i=0; i < *nvar_out; i++) {
                strncpy(record_label, labels[i], 4);
            for(k=0; k<=4; k++) {
              if(record_label[k] == ' ') {
                record_label[k] = '_';
              }
            }
                *(record_label+4) = '\0';
                p = strcat(record_label, min);
                status = ncmpi_put_att_double(ncid, NC_GLOBAL, p, NC_DOUBLE, 1, &varMin[i]);
            if (status < 0){
              printf("Error: ncmpi_put_att_text, varMin\n");
              Driver_abortFlash("Error: ncmpi_put_att_text, varMin\n");
            }
            }

  */
  assert(*pFile_version == FILE_FORMAT_VERSION);
#endif


  /* write out simulation information.  For now, keep these in the global name space */
  file_version = *pFile_version;
  status = ncmpi_put_att_int(ncid, NC_GLOBAL, "file_format_version", NC_INT, 1, &file_version);
  if (status < 0){
    printf("Error: ncmpi_put_att_int, file_format_version\n");
    Driver_abortFlashC("Error: ncmpi_put_att_int, file_format_version\n");
  }

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "setup_call", long_len, setup_call);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, setup_call\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, setup_call\n");
  }

  /*ensure this doesn't change across processors*/
  MPI_Bcast(file_creation_time, strlen(file_creation_time) + 1, MPI_CHAR, MASTER_PE, MPI_COMM_WORLD);

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "file_creation_time", MAX_STRING_LENGTH, file_creation_time);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, file_creation_time\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, file_creation_time\n");
  }

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "flash_version",  MAX_STRING_LENGTH, flash_version);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, flash_version\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, flash_version\n");
  }
  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "build_date", MAX_STRING_LENGTH, build_date); 
  if (status < 0){
    printf("Error: ncmpi_put_att_text, build_date\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, build_date\n");
  }
  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "build_dir", MAX_STRING_LENGTH, build_dir);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, build_dir\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, build_dir\n");
  }
  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "build_machine", MAX_STRING_LENGTH, build_machine);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, build_machine\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, build_machine\n");
  }

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "cflags", long_len, cflags);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, cflags\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, cflags\n");
  }

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "fflags", long_len, fflags);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, fflags\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, fflags\n");
  }


  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "setup_time_stamp", MAX_STRING_LENGTH, setup_time_stamp);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, setup_time_stamp\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, setup_time_stamp\n");
  }

  status = ncmpi_put_att_text(ncid, NC_GLOBAL, "build_time_stamp", MAX_STRING_LENGTH, build_time_stamp);
  if (status < 0){
    printf("Error: ncmpi_put_att_text, build_time_stamp\n");
    Driver_abortFlashC("Error: ncmpi_put_att_text, build_time_stamp\n");
  }









  for(i=0; i < *numIntParms; i++) {
    p = strtok(intParmNames[i], " ");
    status = ncmpi_put_att_int(ncid, 0, p, NC_INT, 1, &intParmValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_int, intParmValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_int, intParmValues\n");
      
    }
  }
  
  
  for(i=0; i < *numRealParms; i++) {
    p = strtok(realParmNames[i], " ");
    status = ncmpi_put_att_double(ncid, 0, p, NC_DOUBLE, 1, &realParmValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_double, realParmValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_double, realParmValues\n");
    }
  }
  

  for(i=0; i < *numStrParms; i++) {
    p = strtok(strParmNames[i], " ");
    status = ncmpi_put_att_text(ncid, 0, p, MAX_STRING_LENGTH, strParmValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_text, strParmValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_text, strParmValues\n");
      
    } 
  }
  
   
  for(i=0; i < *numLogParms; i++) {
    p = strtok(logParmNames[i], " ");
    status = ncmpi_put_att_int(ncid, 0, p, NC_INT, 1, &logToIntParmValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_int, logParmValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_int, logParmValues\n");
      
    }
  }


  for(i=0; i < *numIntScalars; i++) {
    p = strtok(intScalarNames[i], " ");
    status = ncmpi_put_att_int(ncid, 1, intScalarNames[i], NC_INT, 1, &intScalarValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_int, intScalarValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_int, intScalarValues\n");
      
    }
  }

  
  for(i=0; i < *numRealScalars; i++) {
    p = strtok(realScalarNames[i], " ");
    status = ncmpi_put_att_double(ncid, 1, p, NC_DOUBLE, 1, &realScalarValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_double, realScalarValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_double, realScalarValues\n");
    }
  }
  

  for(i=0; i < *numStrScalars; i++) {
    p = strtok(strScalarNames[i], " ");
    status = ncmpi_put_att_text(ncid, 1, p, MAX_STRING_LENGTH, strScalarValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_text, strScalarValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_text, strScalarValues\n");
    } 
  }
  

  for(i=0; i < *numLogScalars; i++) {
    p = strtok(logScalarNames[i], " ");
    status = ncmpi_put_att_int(ncid, 1, p, NC_INT, 1, &logToIntScalarValues[i]);
    if (status < 0){
      printf("Error: ncmpi_put_att_int, logScalarValues\n");
      Driver_abortFlashC("Error: ncmpi_put_att_int, logScalarValues\n");
    }
  }



/* End define mode for standard FLASH I/O.  Do not end define mode
   for experimental I/O because it is closed later in io_writeData.F90 */

#ifndef FLASH_IO_EXPERIMENTAL
  /* end define mode */
  status = ncmpi_enddef(ncid);

  if (status < 0){
    printf("Error: enddef %d\n", status);
#ifdef DEBUG_IO
  if (status != NC_NOERR)
    handle_error(status);
  fflush(stdout);
  fflush(stderr);
#endif
  Driver_abortFlashC("Error: enddef\n");
  }
#endif




}

























