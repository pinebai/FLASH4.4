#include <stdio.h>
#include <string.h>
#include "HYPRE_config.h"
#include "mangle_names.h"
#include "gr_hypre_fortran_message.h"

/* The function prototypes are here because the gr_hypre_fortran_message.h
   header file must be usable by Fortran source files too */
void FTOC(gr_hypre_fortran_message)(char msg[], int *pLen);
void hypre_message(char msg[]);
void hypre_message_test1();
void hypre_message_test2();


/* Called by Fortran */
void FTOC(gr_hypre_fortran_message)(char msg[], int *pLen)
{
  const int len = *pLen;
  int inUseLen;

  if (len >= HYPRE_MSG_LEN) {
    hypre_message(msg);

    /* Add empty characters for Fortran (replace the null character(s)) */
    inUseLen = strlen(msg);
    memset(msg+inUseLen, ' ', len-inUseLen);
  }
}

/* Internal function.  Be careful when calling it directly because the
   error checking is in gr_hypre_fortran_message */
void hypre_message(char msg[])
{
  sprintf(msg,
	  "HYPRE_RELEASE_NAME=%.10s, "
	  "HYPRE_RELEASE_VERSION=%.10s, "
	  "HYPRE_RELEASE_DATE=%.10s",
	  HYPRE_RELEASE_NAME,
	  HYPRE_RELEASE_VERSION,
	  HYPRE_RELEASE_DATE);
}



/*
  To build the hypre message unit test with gcc:

  gcc gr_hypre_fortran_message.c     \
  -I /opt/hypre/gnu/2.7.0b/include   \
  -I ../../../flashUtilities/general \
  -DHYPRE_MESSAGE_TEST -o unit_test

  (Note that there is no need to link against HYPRE because the unit
  test only use the HYPRE macro definitions.)
*/

#ifdef HYPRE_MESSAGE_TEST
int main()
{
  printf("In the test we pretty-print HYPRE information:\n");
  hypre_message_test1();
  hypre_message_test2();
  printf("SUCCESS!\n");
  return 0;
}

void hypre_message_test1()
{
  char msg[HYPRE_MSG_LEN];
  hypre_message(msg);
  printf("Test1: %s\n", msg);
}

void hypre_message_test2()
{
  char msg[HYPRE_MSG_LEN];
  int len = HYPRE_MSG_LEN;
  FTOC(gr_hypre_fortran_message)(msg, &len);
  msg[HYPRE_MSG_LEN-1] = '\0';
  printf("Test2: %s\n", msg);
}
#endif
