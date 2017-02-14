#ifndef _OPTIONS_H
#define _OPTIONS_H

#include "constants.h"

typedef struct options_t{
  int splitNum; /*number of split files to make one whole file*/
  int zeroPad;  /*the total size of digits for the splitnum*/
  char basenm[MAX_STRING_LENGTH]; /*basename for the file*/
  char format[10]; /*file format (hdf5, pnetcdf, etc.)*/
  char type[20]; /*plot file, particle file or checkpoint for now*/
  char filenum[10]; /*the file number of the overall file*/
  /*char* filename; /*target filename (may not be used in the end*/
  char filename[STRINGSIZE] 
}options_t;

int parse_cmdline(options_t *opts, int argc, char **argv);

#endif /* _OPTIONS_H */
