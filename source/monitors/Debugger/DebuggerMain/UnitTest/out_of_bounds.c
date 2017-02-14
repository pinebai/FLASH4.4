#include <stdlib.h>
#include <stdio.h>
#include "mangle_names.h"
#include "dbg_heap_check.h"

int main(int argc, char *argv[])
{
  double *x;
  const size_t size=10;
  size_t i;

  x = malloc(size*sizeof(*x));
  if (NULL != x) {
    for (i=0; i<size+1; ++i) {
      x[i] = 12.9;
      printf("x[%zu]=%lf\n", i, x[i]);
    }
    FTOC(dbg_heap_check)();
    free(x);
  }

  printf("The memory corruption was not detected. FAIL!\n");
  return 0;
}
