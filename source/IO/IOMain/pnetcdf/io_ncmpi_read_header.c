/* This file contains the functions that read the data from the pnetcdf file
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
#include <pnetcdf.h>
#include "mangle_names.h"
#include <mpi.h>
#include <assert.h>
#include "Flash.h"
#include "constants.h"


int Driver_abortFlashC(char* message);


/* reads runtime parameter and scalar data */

void FTOC(io_ncmpi_read_header)(int* MyPE,
                          int* file_identifier,    /* file handle */
                          int* nvar_out,
                          char unk_labels[][4],
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
                          int logToIntScalarValues[])
   
{




  int ncid, status, i, k;
  char record_label[5];
  char temp_name[MAX_STRING_LENGTH];

  int varid[40];
  int id = 0;
  int rp_id = 0, s_id = 1;
  int nruntime_parms;
  int nscalars;
  int temp_int;
  double temp_real;
  nc_type xytpep; /*nc_type returned by inq_atttype */


  ncid = *file_identifier;
  *numRealParms = 0;
  *numIntParms = 0;
  *numStrParms = 0;
  *numLogParms = 0;
  *numRealScalars = 0;
  *numIntScalars = 0;
  *numStrScalars = 0;
  *numLogScalars = 0;


  /*attempt to get the runtime parameters*/

  status = ncmpi_inq_varnatts(ncid, rp_id, &nruntime_parms);


  for(i=0; i< nruntime_parms; i++){
    status = ncmpi_inq_attname(ncid, rp_id, i, temp_name);
    status = ncmpi_inq_atttype(ncid, rp_id, temp_name, &xytpep);

    if (xytpep == NC_INT){
      /* a bit convoluted but necessary because we
       store logical values as ints.  This if 
       statement just checks verifies we are reading in 
       either ints or the logical values stored as ints */

      if (*numRealParms == 0){
      status = ncmpi_get_att_int(ncid, rp_id, temp_name, &temp_int);
      intParmValues[*numIntParms] = temp_int;
      strcpy(intParmNames[*numIntParms], temp_name);
      (*numIntParms) ++;
      }else{
      status = ncmpi_get_att_int(ncid, rp_id, temp_name, &temp_int);
      logToIntParmValues[*numLogParms] = temp_int;
      strcpy(logParmNames[*numLogParms], temp_name);
      (*numLogParms) ++;
      }

    }else if (xytpep == NC_DOUBLE){
      status = ncmpi_get_att_double(ncid, rp_id, temp_name, &temp_real);
      if(status < 0){
      Driver_abortFlashC("error getting real att\n");
      }
      realParmValues[*numRealParms] = temp_real;
      strcpy(realParmNames[*numRealParms], temp_name);
      (*numRealParms) ++;

    }else if (xytpep == NC_CHAR){
      status = ncmpi_get_att_text(ncid, rp_id, temp_name, strParmValues[*numStrParms]);
      strcpy(strParmNames[*numStrParms], temp_name);
      (*numStrParms) ++;
    }else{
      Driver_abortFlashC("Don't know what kind of attribute you are attempting to read in\n");
    }
    

  }



  /*attempt to get the scalars*/

  status = ncmpi_inq_varnatts(ncid, s_id, &nscalars);


  for(i=0; i< nscalars; i++){
    status = ncmpi_inq_attname(ncid, s_id, i, temp_name);
    status = ncmpi_inq_atttype(ncid, s_id, temp_name, &xytpep);

    if (xytpep == NC_INT){
      /* a bit convoluted but necessary because we
       store logical values as ints.  This if 
       statement just checks verifies we are reading in 
       either ints or the logical values stored as ints */

      if (*numRealScalars == 0){
      status = ncmpi_get_att_int(ncid, s_id, temp_name, &temp_int);
      strcpy(intScalarNames[*numIntScalars], temp_name);
      intScalarValues[*numIntScalars] = temp_int;
      (*numIntScalars) ++;
      }else{
      status = ncmpi_get_att_int(ncid, s_id, temp_name, &temp_int);
      strncpy(logScalarNames[*numLogScalars], temp_name, MAX_STRING_LENGTH);
      logToIntScalarValues[*numLogScalars] = temp_int;
      (*numLogScalars) ++;
      }

    }else if (xytpep == NC_DOUBLE){
      status = ncmpi_get_att_double(ncid, s_id, temp_name, &temp_real);
      if(status < 0){
      Driver_abortFlashC("error getting real att\n");
      }
      strncpy(realScalarNames[*numRealScalars], temp_name, MAX_STRING_LENGTH);
      realScalarValues[*numRealScalars] = temp_real;
      (*numRealScalars) ++;

    }else if (xytpep == NC_CHAR){
      status = ncmpi_get_att_text(ncid, s_id, temp_name, strScalarValues[*numStrScalars]);
      strncpy(strScalarNames[*numStrScalars], temp_name, MAX_STRING_LENGTH);
      (*numStrScalars) ++;
    }else{
      Driver_abortFlashC("Don't know what kind of attribute you are attempting to read in\n");
    }
  }



#ifdef DEBUG_IO
  if (status != NC_NOERR)  handle_error(status);
#endif
  
}





