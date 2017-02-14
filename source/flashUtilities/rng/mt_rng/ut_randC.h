#ifdef UNIT_TEST_MODE
#include <stdio.h>
#endif

#include <stdlib.h>
#include "mangle_names.h"

void FTOC(ut_rand_seed)(int *getseed, int *putSeed);
void FTOC(ut_rand)(double *x);
void FTOC(ut_rand_init)();

