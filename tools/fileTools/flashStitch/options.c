#include <unistd.h>
#include <stdlib.h>
#include "options.h"

int parse_cmdline(options_t *opts, int argc, char **argv){
  
  char *file1, *file2, *pwd;
  int i,j,pwdLen, filenameLen;
  int argError=0, c;
  int fileNum;

  /* initialize opts */
  opts->splitNum = 0;
  opts->zeroPad = 4; 
  strncpy(opts->format, "hdf5", 5);
  strncpy(opts->type, "chk", 5);

  
  /* Argument parsing */
  while ((c = getopt(argc, argv, "s:z:n:f:t:i:")) != -1) {
    switch(c) {
    case 's':
      opts->splitNum = atoi(optarg);
      break;
    case 'z':
      opts->zeroPad = atoi(optarg);
      break;
    case 'n':
      strncpy(opts->basenm, optarg, MAX_STRING_LENGTH);
      break;
    case 'f':
      strncpy(opts->format, optarg, 20);
      break;
    case 't':
      strncpy(opts->type, optarg, 10);
      break;
    case 'i':
      strncpy(opts->filenum, optarg, 10);
      break;
    case '?':
      /* do nothing - this is a change to prevent sfocu from
         crashing when it tries to interpret options meant for
         mpirun or poe. I have left the original 'return 1' in
         this commented-out block
      return 1;
      Note form PR: I have also left this in for future mpi use.
      */
      break;
    default: 
      printf("Unrecognized option.  Exiting...\n");
      return -1;
      break;
    }
  }

  /*final three things should be required arguments*/
  if(argc < 3)
    /*we have a problem here*/
    return -1;
  
  strncpy(opts->basenm, argv[argc-3], MAX_STRING_LENGTH);
  opts->splitNum = atoi(argv[argc-2]);
  strncpy(opts->filenum, argv[argc-1], 10);
  

  /*filename is last argument*/
  sprintf(opts->filename, "%s_%s_%s_%s", opts->basenm, opts->format,
          opts->type, opts->filenum);

  return 0;
}
