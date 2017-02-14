#ifndef UT_SYSMEMPROC_H
#define UT_SYSMEMPROC_H

#include <stdio.h>
#include "ut_sysMemCTypes.h"

/* 04/17/2012 - Avoid using /proc/self/stat on IBM machines.  This
   is the way FLASH has done things for a long time - see io_memory_usage.c
   from <= FLASH4-beta.  "IBM" is a FLASH Makefile.h defined macro and
   "__bg__" is a Blue Gene compiler defined macro. */
#if defined(__linux__) && !(defined(IBM) || defined(__bg__))
# define FLASH_SUPPORT_PROC
#endif

#ifdef FLASH_SUPPORT_PROC
# define PROC_STATS 2
# include <unistd.h>
#else
# define PROC_STATS 0
#endif

int ut_sysMemProcNumStats(int verbosity);
void ut_sysMemProc(meminfo_t *meminfo, int meminfoSize, int verbosity);

#endif
