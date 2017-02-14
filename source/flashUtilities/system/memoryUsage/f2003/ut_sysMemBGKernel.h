#ifndef UT_SYSMEMBGKERNEL_H
#define UT_SYSMEMBGKERNEL_H

/*
   It is possible that the header file containing the function prototype
   for Kernel_GetMemorySize is not in our include path.  I have noticed this
   on BG/P but not on BG/Q.

   If we are using the default MPI wrapper scripts and the BG/P sys-admins
   have kept the software up to date then we just need to add the following
   include path to a CFLAGS macro in Makefile.h
   > -I/bgsys/drivers/ppcfloor/arch/include

   The /bgsys/drivers path is a soft link
   > ls -l /bgsys/drivers/
   > ppcfloor -> /bgsys/drivers/V1R4M2_200_2010-100508P/ppc

   We can check whether our MPI wrapper script was built using the same
   driver version by typing
   > mpixlc_r -show

   This will return something like
   > ... -I/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/default/include ...
   which agrees with the ppcfloor soft link
*/

#include "ut_sysMemCTypes.h"

#ifdef __bg__
# define FLASH_SUPPORT_BGKERNEL
#endif

#ifdef FLASH_SUPPORT_BGKERNEL
/* BG/Q has 64-bit processors.  BG/L and BG/P have 32-bit processors. */
# ifdef __PPC64__
#  include <spi/include/kernel/memory.h>
typedef uint64_t BG_memsize_int;
# else
#  include <spi/kernel_interface.h>
typedef uint32_t BG_memsize_int;
# endif
#endif

int ut_sysMemBGKernelNumStats(int verbosity);
void ut_sysMemBGKernel(meminfo_t *meminfo, int meminfoSize, int verbosity);

#endif
