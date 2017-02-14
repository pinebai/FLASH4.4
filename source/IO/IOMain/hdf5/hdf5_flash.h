/* general header file for the HDF 5 IO in FLASH */


#ifndef _HDF5_FLASH_H
#define _HDF5_FLASH_H

/* pull in some basic FLASH information */


/* define an integer file format version number, that is stored
   in the output files.  This way, people can check this number
   before reading, and compare against a published format to know
   what is stored in the file.  In theory, this number should be
   incremented anytime a change is made to the file format

   4 -- extrema attributes added

   5 -- redshift included (incremented to agree with serial version)

   6 -- added the Module data to the attributes of "/"
   
   7 -- made build info attributes instead of datasets 

   8 -- flash3 IO.  Differences handling scalars, also
        no longer has sim_params structure holding time, dt, nxb, etc
	-- These values are now stored in the scalars list and parameters list

   9 -- Fields that were formerly of size NDIM are now MDIM, so that data is
        not lost when Flash uses Paramesh in the 2.5d mode.

*/

#include "constants.h"


#define FILE_FORMAT_VERSION 9

#define LIST_STRING_SIZE MAX_STRING_LENGTH

typedef struct int_list_t {
  char name[LIST_STRING_SIZE];
  int value;
} int_list_t;


typedef struct real_list_t {
  char name[LIST_STRING_SIZE];
  double value;
} real_list_t;


typedef struct str_list_t {
  char value[LIST_STRING_SIZE];
  char name[LIST_STRING_SIZE];
} str_list_t;


typedef struct log_list_t {
  int value;
  char name[LIST_STRING_SIZE];
} log_list_t;



typedef struct sim_info_t {
  int file_format_version;
  char setup_call[400];
  char file_creation_time[MAX_STRING_LENGTH];
  char flash_version[MAX_STRING_LENGTH];
  char build_date[MAX_STRING_LENGTH];
  char build_dir[MAX_STRING_LENGTH];
  char build_machine[MAX_STRING_LENGTH];
  char cflags[400];
  char fflags[400];
  char setup_time_stamp[MAX_STRING_LENGTH];
  char build_time_stamp[MAX_STRING_LENGTH];
} sim_info_t;





/* set the dimension and grid variables -- the variable N_DIM is set 
   in the compile line */


/* 3-d problem */
#if N_DIM == 3 

#define NDIM  3


#define NGID 15

#define k2d 1
#define k3d 1


/* 2-d problem */
#elif N_DIM == 2

#define NDIM  2


#define NGID 9

#define k2d 1
#define k3d 0


/* 1-d problem */
#else

#define NDIM 1

#define NGID 5

#define k2d 0
#define k3d 0

#endif


/*mflags hard coded as 1 for bflags data struct in PM3 */
#define MFLAGS 1

/*This will let us know if we are using collective mode for our HDF5
  write calls (as opposed to independent mode, which is default).
  This only really makes sense when HDF5 is run in parallel mode.
  Allowing flexibility for use with other modes (compression, chunking, etc.)
*/
extern int HDF5_MODE;
#define INDEPENDENT 0
#define COLLECTIVE 1

#define MAX_ARR_DIMS 5

#endif

