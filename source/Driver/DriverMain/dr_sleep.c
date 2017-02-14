#include <unistd.h>
#include <stdio.h>
#include "mangle_names.h"

void FTOC(dr_sleep)(int *sleepInt)
{
  unsigned int sleepUint, rtn;

  if (*sleepInt >= 0)
  {
    sleepUint = (unsigned int) *sleepInt;
    rtn = sleep(sleepUint);
  
#ifdef DEBUG_ALL
    if (rtn != 0)
      fprintf(stderr, "[dr_sleep]: Sleep interupted..."
              "time remaining: %u seconds\n", rtn);
#endif
  }
}
