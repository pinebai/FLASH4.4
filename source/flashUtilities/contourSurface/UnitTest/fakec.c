
#include <stdio.h>
#include <stdlib.h>

int Driver_abortFlashC(char* message){
  fprintf(stderr, message);
  exit(1);

  return 0;
}

