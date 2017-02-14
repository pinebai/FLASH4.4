#ifndef FLASHSTITCHER_H
#define FLASHSTITCHER_H

#include "hdf5.h"
#include "options.h"

/*
#define MAX_STRING_LENGTH 80
#define STRINGSIZE 256
#define UNK_NAME_LEN 4
*/

/*custom types*/
typedef struct int_list_t{
  char name[MAX_STRING_LENGTH];
  int value;
} int_list_t;

typedef struct real_list_t{
  char name[MAX_STRING_LENGTH];
  double value;
} real_list_t;
  
typedef struct log_list_t{
  char name[MAX_STRING_LENGTH];
  int value; /*use int_list_t?*/
} log_list_t;

typedef struct str_list_t{
  char name[MAX_STRING_LENGTH];
  char value[MAX_STRING_LENGTH];
} str_list_t;

/*metadata type*/
typedef struct sim_info_t{
  int  file_format_version;
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
}sim_info_t;

/*function definitions*/
int flashStitcher(char* filename, int splitNum, int zeroPad, options_t *opts);
int readFlashInfo(hid_t *file,sim_info_t *siminfo, char **unkNames, int *numUnkNames);
int writeFlashInfo(hid_t * file, sim_info_t *siminfo, char *unkNames, int *numUnkNames);
void printSimInfo(sim_info_t *si);
void copy_bflags(hid_t infile, hid_t outfile, int currentOffset, 
                 int globalNumBlocks, int *nextOffset);
void copy_blocksize(hid_t infile, hid_t outfile, int currentOffset, 
                    int globalNumBlocks);
void copy_coordinates(hid_t infile, hid_t outfile, int currentOffset, 
                    int globalNumBlocks);
void copy_gid(hid_t infile, hid_t outfile, int currentOffset, 
                    int globalNumBlocks);
void copy_bndbox(hid_t infile, hid_t outfile, int currentOffset, 
                    int globalNumBlocks);
void copy_1d_dataset(hid_t infile, hid_t outfile, int currentOffset, 
                     int globalNumBlocks, hid_t type, char *name, 
                     int firstCall);
void copy_unkvars(hid_t infile, hid_t outfile, int currentOffset, 
                  int globalNumBlocks, char *unkNames, int numUnkNames,
                  int firstCall);
#endif
