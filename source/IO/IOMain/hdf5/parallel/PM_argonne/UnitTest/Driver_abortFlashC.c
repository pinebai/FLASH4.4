#include <stdio.h>
#include <stdlib.h>

/* Quick implementation for serial unit tests */

int Driver_abortFlashC(const char message[])
{
  fprintf(stderr, "%s\n", message);
  exit(EXIT_FAILURE);
  return EXIT_FAILURE;
} 
