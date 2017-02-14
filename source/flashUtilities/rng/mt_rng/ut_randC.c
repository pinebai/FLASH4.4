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

#include "ut_randC.h"
#include "mt19937ar.h"

void FTOC(ut_rand_init)()
{
  init_genrand(5489UL);
  mti=N+1;
}

void FTOC(ut_get_size)(int *t_size)
{
  *t_size=(N+2)*sizeof(unsigned long);
}
  
void FTOC(ut_get_seed)(int *vec)
{
  unsigned long *temp_vec;
  int i;
  temp_vec=(unsigned long *)vec;
  i=M;
  temp_vec[0]=(unsigned long)i;
  temp_vec[1]=(unsigned long)mti;
  for(i=0;i<N;i++)
    {
      temp_vec[i+2]=mt[i];
    }
}

void FTOC(ut_put_seed)(int *vec)
{
  unsigned long *temp_vec;
  int i,j;

  temp_vec=(unsigned long *)vec;
  j=(int)temp_vec[0];
  printf("the value in j %d \n",j);
  if(j==M)
    {
      mti=(int)temp_vec[1];
      printf("the value in mti %d \n",j);
      if(mti<N+1) 
	{
	  for(i=0;i<N;i++)
	    {
	      mt[i]=temp_vec[i+2];
	    }
	}
    }
  else
    {
      mti=N+1;
    }
}

/* generates a random number on [0,1) with 53-bit resolution*/
void FTOC(ut_rand)(double *x)
{
  *x=genrand_real1();
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


