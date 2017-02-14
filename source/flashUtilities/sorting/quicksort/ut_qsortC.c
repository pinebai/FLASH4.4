/*  
    No function overloading in C.  I could write it in C++ with 
    function templates but then it is hard to know mangled function name, 
    and we could not call easily from Fortran

    http://www.math.utah.edu/software/c-with-fortran.html
    For mixed-language programming in C/C++ and Fortran, only the
    minimal intersection of their many data types can be relied on:

    * int == INTEGER,
    * float == REAL,
    * double == DOUBLE PRECISION.
    
    No other data types can be expected to be exchanged without serious
    compromise of portability.
*/

#include "ut_qsortC.h"
  

/* Integer sorting wrappers */
void FTOC(ut_qsort_int_asc)(int inOutArray[], 
			    const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(int), compareIntAsc);
}

void FTOC(ut_qsort_int_desc)(int inOutArray[], 
			     const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(int), compareIntDesc);
}


/* Float sorting wrappers */
void FTOC(ut_qsort_float_asc)(float inOutArray[], 
			      const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(float), compareFloatAsc);
}

void FTOC(ut_qsort_float_desc)(float inOutArray[], 
			       const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(float), compareFloatDesc);
}


/* Double sorting wrappers */
void FTOC(ut_qsort_double_asc)(double inOutArray[], 
			       const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(double), compareDoubleAsc);
}

void FTOC(ut_qsort_double_desc)(double inOutArray[], 
				const int * const pArrLength)
{
  assert(*pArrLength > 0);
  qsort(inOutArray, (size_t) *pArrLength, sizeof(double), compareDoubleDesc);
}



/* Comparison functions for qsort stdlib.h function */
int compareIntAsc(const void *a, const void *b)
{
  return (*(int*)a - *(int*)b);
}

int compareIntDesc(const void *a, const void *b)
{
  return (*(int*)b - *(int*)a);
}

int compareFloatAsc(const void *a, const void *b)
{
  return (*(float*)a - *(float*)b);
}

int compareFloatDesc(const void *a, const void *b)
{
  return (*(float*)b - *(float*)a);
}

int compareDoubleAsc(const void *a, const void *b)
{
  return (*(double*)a - *(double*)b);
}

int compareDoubleDesc(const void *a, const void *b)
{
  return (*(double*)b - *(double*)a);
}



#ifdef UNIT_TEST_MODE
int main()
{
#define SIZE 5
  const int sizeArr = SIZE;
  int intArray[SIZE] = {4,13,7,82,1};
  float floatArray[SIZE] = {4.1,13.1,7.1,82.1,1.1};
  double doubleArray[SIZE] = {4.1,13.1,7.1,82.1,1.1};
  int i;

  printf(" **** ORIGINAL INTEGER ARRAY ****\n");
  for (i=0; i<sizeArr; ++i) {
    printf("intArray[%d]=%d\n", i, intArray[i]);
  }
  
  FTOC(ut_qsort_int_asc)(intArray, &sizeArr);

  printf("\nASCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("intArray[%d]=%d\n", i, intArray[i]);
  }
  
  FTOC(ut_qsort_int_desc)(intArray, &sizeArr);

  printf("\nDESCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("intArray[%d]=%d\n", i, intArray[i]);
  }



  printf("\n\n\n **** ORIGINAL FLOAT ARRAY ****\n");
  for (i=0; i<sizeArr; ++i) {
    printf("floatArray[%d]=%8.1f\n", i, floatArray[i]);
  }
  
  FTOC(ut_qsort_float_asc)(floatArray, &sizeArr);

  printf("\nASCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("floatArray[%d]=%8.1f\n", i, floatArray[i]);
  }
  
  FTOC(ut_qsort_float_desc)(floatArray, &sizeArr);

  printf("\nDESCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("floatArray[%d]=%8.1f\n", i, floatArray[i]);
  }



  printf("\n\n\n **** ORIGINAL DOUBLE ARRAY ****\n");
  for (i=0; i<sizeArr; ++i) {
    printf("doubleArray[%d]=%8.1lf\n", i, doubleArray[i]);
  }
  
  FTOC(ut_qsort_double_asc)(doubleArray, &sizeArr);

  printf("\nASCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("doubleArray[%d]=%8.1lf\n", i, doubleArray[i]);
  }
  
  FTOC(ut_qsort_double_descC)(doubleArray, &sizeArr);

  printf("\nDESCENDING ORDER\n");
  for (i=0; i<sizeArr; ++i) {
    printf("doubleArray[%d]=%8.1f\n", i, doubleArray[i]);
  }
}
#endif
