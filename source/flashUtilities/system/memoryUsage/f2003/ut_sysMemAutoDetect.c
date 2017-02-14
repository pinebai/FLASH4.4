#include "ut_sysMemAutoDetect.h"

int ut_sysMemAutoDetect(void)
{
  int memorySampler = 0;

  /* On linux machines we take memory usage information from /proc/self/stat. */
#ifdef __linux__
  memorySampler = UT_SYSMEM_PROC;
#endif

  /* On IBM machines we prefer to use rusage instead of /proc/self/stat.
     IBM is a FLASH-defined macro. */
#ifdef IBM
  memorySampler = UT_SYSMEM_RUSAGE;
#endif

#ifdef __bg__
# ifdef __PPC64__
  /* On BG/Q we use bgkernel because neither /proc/self/stat nor rusage work. */
  memorySampler = UT_SYSMEM_BGKERNEL;
# else
  /* On BG/P we use both rusage and bgkernel. */
  memorySampler = UT_SYSMEM_RUSAGE | UT_SYSMEM_BGKERNEL;
# endif
#endif

  return memorySampler;
}
