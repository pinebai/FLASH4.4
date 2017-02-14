#ifndef UT_SYSMEMMALLINFO_H
#define UT_SYSMEMMALLINFO_H

#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include "ut_sysMemCTypes.h"

#if defined (__GLIBC__) || defined (__GNU_LIBRARY__)
# define FLASH_SUPPORT_MALLINFO
#endif

int ut_sysMemMallinfoNumStats(int verbosity);
void ut_sysMemMallinfo(meminfo_t *meminfo, int meminfoSize, int verbosity);

#endif
