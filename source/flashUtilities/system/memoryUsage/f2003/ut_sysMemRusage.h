#ifndef UT_SYSMEMRUSAGE_H
#define UT_SYSMEMRUSAGE_H

#include <stdio.h>
#include "ut_sysMemCTypes.h"

/* 04/17/2012 - Avoid using getrusage on BG/Q (__bg__ && __PPC64__)
   because it is returning garbage.  This may just be a short-term
   problem which is happening because BG/Q is currently in alpha state */
#if defined(__unix__) && !(defined(__bg__) && defined(__PPC64__))
# define FLASH_SUPPORT_RUSAGE
#endif

#ifdef FLASH_SUPPORT_RUSAGE
# define RUSAGE_STATS 1
# include <sys/time.h>
# include <sys/resource.h>
#else
# define RUSAGE_STATS 0
#endif

int ut_sysMemRusageNumStats(int verbosity);
void ut_sysMemRusage(meminfo_t *meminfo, int meminfoSize, int verbosity);

#endif
