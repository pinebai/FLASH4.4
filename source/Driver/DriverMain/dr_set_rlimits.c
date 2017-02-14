#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include "mangle_names.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "constants.h"

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

void FTOC(dr_set_rlimits)(int* myPE) {

  int retval;
  int before_cur;
  struct rlimit lim;
  mode_t err;

  mode_t mask=07; /* this bit is to make sure that the files are created with permissions 
                     we want and not what the system wants */
  err=umask(mask);

  retval = getrlimit(RLIMIT_STACK, &lim);
  before_cur = lim.rlim_cur;

  lim.rlim_cur = RLIM_INFINITY;
  lim.rlim_max = RLIM_INFINITY;

  retval = setrlimit(RLIMIT_STACK, &lim);

  retval = getrlimit(RLIMIT_STACK, &lim);

#ifdef DEBUG_DRIVER
  if (before_cur != lim.rlim_cur) {
    if (*myPE == MASTER_PE)
      printf (" dr_set_rlimits: changed current stack limit from: %d to %d\n", before_cur, lim.rlim_cur);
  }
#endif
}
