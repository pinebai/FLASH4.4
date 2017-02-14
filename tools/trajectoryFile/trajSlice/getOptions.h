#ifndef __GET_OPTIONS_H
#define __GET_OPTIONS_H

#include <getopt.h>

typedef struct arguments_t{
  char filenameIn[256];
  char filenameOut[256];
  int  numToRead;
  int  start;
  int  stride;
  int  useThreshold;
  int  findAllThreshold;
  char thresholdVarName[24];

} arguments_t;

//arguments_t arguments;

static struct option long_options[] = 
  {
    {0,0,0,0}
  };


#define ARGUMENT_ERROR -1
#define OPTION_ERROR -2

int getOptions(int* argc, char** argv, arguments_t * arguments);

#endif
