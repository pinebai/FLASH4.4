/*get command line options and process command line arguments*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "getOptions.h"

int getOptions(int* argc, char **argv, arguments_t* arguments){


  int option_index = 0;
  int c = -1;
  //char* option_index;
  
  if(*argc < 5){ /*insufficient arguments have been passed*/
    printf("Insufficient arguments\n");
    return ARGUMENT_ERROR;
  }
   
  strncpy(arguments->filenameOut, argv[*argc - 1], 256);
  strncpy(arguments->filenameIn, argv[*argc - 2], 256);

  arguments->stride  = atoi(argv[*argc-3]);
  arguments->start  = atoi(argv[*argc-4]);
  arguments->numToRead = atoi(argv[*argc-5]);

  arguments->useThreshold = 0;
  arguments->findAllThreshold = 0;
  

  while(1){
    
    if((c = getopt_long(*argc, argv, "ft", long_options, &option_index)) == -1)
      break;
    switch (c){
    case 't': 
      arguments->useThreshold = 1;
      break;
    case 'f':
      arguments->findAllThreshold = 1;
      break;
     case '?':
      printf ("unknown option detected\n");
      break;
    default:
      printf("should never be here!");
      return OPTION_ERROR;
    }

  }
  
  strncpy(arguments->thresholdVarName, "flam", 24);

  // printf("argument values:\n pad %d split %d flush %d\n %s %d %d\n",
  //     arguments->padDigits, arguments->splitNumber,
  //     arguments->particlesPerFlush, arguments->basename,
  //     arguments->start_num, arguments->end_num);

  return 0;
}
