#if 0
This file is included by Fortran and C++ files.

It specifies integer values for true and false because
of portability issues with interoperability of boolean types.

A logical(c_bool) is interoperable with _Bool which is part
of the C99 standard but not part of the C++ standard.

Our strategy for interoperability with C++ will be to use
integer(c_int) with the convention that true is 1 and
false is 0

The defines are named FLASH_TRUE and FLASH_FALSE instead of
TRUE and FALSE because application code may have .TRUE.
and this would become .1. after preprocessing.
#endif

#define FLASH_TRUE 1
#define FLASH_FALSE 0
