#if 0
This header file provides HYPRE_MSG_LEN to C and Fortran source files
HYPRE_MSG_LEN is the minimum size of the character buffer which you
can pass to gr_hypre_fortran_message
#endif

#define NAME_LEN (18 + 1 + 10 + 2)
#define VERSION_LEN (21 + 1 + 10 + 2)
#define DATE_LEN (18 + 1 + 10)
#define HYPRE_MSG_LEN (NAME_LEN + VERSION_LEN + DATE_LEN + 1)
