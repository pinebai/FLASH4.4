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
  
  if(*argc < 4){ /*insufficient arguments have been passed*/
    printf("Insufficient arguments\n");
    return ARGUMENT_ERROR;
  }
   
  strncpy(arguments->basename, argv[*argc - 3], 256);
  
  arguments->start_num  = atoi(argv[*argc-2]);
  arguments->end_num  = atoi(argv[*argc-1]);

  arguments->padDigits = 4;
  arguments->splitNumber = 32;
  arguments->timestepsPerFlush = 100;
  
/*  while(1){
    
    if((c = getopt_long(*argc, argv, "p:s:f:", long_options, &option_index)) == -1)
      break;
    switch (c){
    case 'p': //number padding
      arguments->padDigits = atoi(optarg);
      break;
    case 's': //split the file
      arguments->splitNumber = atoi(optarg);
      break;
    case 'f': //number of particles to flush to file
      arguments->timestepsPerFlush = atoi(optarg);
      break;
    case '?':
      printf ("unknown option detected\n");
      break;
    default:
      printf("should never be here!");
      return OPTION_ERROR;
    }

  }
*/
//  arguments->splitNumber = 10;
  // printf("argument values:\n pad %d split %d flush %d\n %s %d %d\n",
  //     arguments->padDigits, arguments->splitNumber,
  //     arguments->particlesPerFlush, arguments->basename,
  //     arguments->start_num, arguments->end_num);

  return 0;
}
