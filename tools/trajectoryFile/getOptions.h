#ifndef __GET_OPTIONS_H
#define __GET_OPTIONS_H

#include <getopt.h>
typedef struct arguments_t{
  char basename[256];
  int start_num;
  int end_num;
  int padDigits;
  int splitNumber;
  int timestepsPerFlush;

} arguments_t;

//arguments_t arguments;

static struct option long_options[] = 
  {
    {"pad",optional_argument,NULL,'p'},
    {"split",required_argument,NULL,'s'},
    {"flush",required_argument,NULL,'f'},
    {0,0,0,0}
  };


#define ARGUMENT_ERROR -1
#define OPTION_ERROR -2

int getOptions(int* argc, char** argv, arguments_t * arguments);

#endif
