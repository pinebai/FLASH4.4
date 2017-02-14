#include "ut_sysMemBGKernel.h"

#ifdef FLASH_SUPPORT_BGKERNEL
typedef struct opmap_t
{
  int verbosity;
  enum KERNEL_MEMSIZETYPE op;
  const char *description;
} opmap_t;

/* /bgsys/drivers/ppcfloor/spi/include/kernel/memory.h on BG/Q and
   /bgsys/drivers/ppcfloor/arch/include/spi/kernel_interface.h on BG/P */

/* All measurements from Kernel_GetMemorySize are in bytes */
static const struct opmap_t opmap[] =
{
  {
    0,
    KERNEL_MEMSIZE_HEAP,
    "bg heap use    (MB):"
  },
  {
    0,
    KERNEL_MEMSIZE_HEAPMAX,
    "bg heap peak   (MB):"
  },
  {
    0,
    KERNEL_MEMSIZE_HEAPAVAIL,
    "bg heap free   (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_STACK,
    "bg stack use   (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_STACKAVAIL,
    "bg stack free  (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_ESTHEAPAVAIL,
    "~bg heap free  (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_ESTSTACK,
    "~bg stack use  (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_ESTSTACKAVAIL,
    "~bg stack free (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_PERSIST,
    "bg persistent  (MB):"
  },
# ifdef __PPC64__
  {
    1,
    KERNEL_MEMSIZE_MMAP,
    "bg mmap        (MB):"
  },
# endif
  {
    1,
    KERNEL_MEMSIZE_GUARD,
    "bg guard page  (MB):"
  },
  {
    1,
    KERNEL_MEMSIZE_SHARED,
    "bg shared      (MB):"
  },
};
# define BGKERNEL_STATS sizeof opmap / sizeof opmap[0]
#endif


int ut_sysMemBGKernelNumStats(int verbosity)
{
  int numMeasurements, i;
  numMeasurements = 0;
#ifdef FLASH_SUPPORT_BGKERNEL
  for (i=0; i<BGKERNEL_STATS; ++i) {
    if (opmap[i].verbosity <= verbosity) {
      ++numMeasurements;
    }
  }
#endif
  return numMeasurements;
}

void ut_sysMemBGKernel(meminfo_t *meminfo, int meminfoSize, int verbosity)
{
#ifdef FLASH_SUPPORT_BGKERNEL
  BG_memsize_int memory_size = 0; /* See typedef in ut_sysMemBGKernel.h */
  const double convertToMB = 1.0 / (1024.0 * 1024.0);
  int numMeasurements, i, j;

  numMeasurements = 0;
  for (i=0; i<BGKERNEL_STATS; ++i) {
    if (opmap[i].verbosity <= verbosity) {
      ++numMeasurements;
      if (numMeasurements <= meminfoSize) {
	j = numMeasurements - 1;
	Kernel_GetMemorySize(opmap[i].op, &memory_size);
	meminfo[j].measurement = memory_size * convertToMB;
	meminfo[j].description = opmap[i].description;
      }
    }
  }
#endif
}
