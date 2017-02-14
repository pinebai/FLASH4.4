#ifndef UT_QSORTC_H
#define UT_QSORTC_H

#ifdef UNIT_TEST_MODE
#include <stdio.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include "mangle_names.h"

void FTOC(ut_qsort_int_asc)(int inOutArray[], 
			    const int * const pArrLength);
void FTOC(ut_qsort_int_desc)(int inOutArray[], 
			     const int * const pArrLength);

void FTOC(ut_qsort_float_asc)(float inOutArray[], 
			      const int * const pArrLength);
void FTOC(ut_qsort_float_desc)(float inOutArray[], 
			       const int * const pArrLength);

void FTOC(ut_qsort_double_asc)(double inOutArray[],
			       const int * const pArrLength);
void FTOC(ut_qsort_double_desc)(double inOutArray[], 
				const int * const pArrLength);

int compareIntAsc(const void *a, const void *b);
int compareIntDesc(const void *a, const void *b);

int compareFloatAsc(const void *a, const void *b);
int compareFloatDesc(const void *a, const void *b);

int compareDoubleAsc(const void *a, const void *b);
int compareDoubleDesc(const void *a, const void *b);

#endif
