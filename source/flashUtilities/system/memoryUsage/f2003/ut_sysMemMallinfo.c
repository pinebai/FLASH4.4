#include "ut_sysMemMallinfo.h"

#ifdef FLASH_SUPPORT_MALLINFO
typedef struct mallinfo_info_t
{
  int verbosity;
  const char *description;
} mallinfo_info_t;

static const struct mallinfo_info_t mallinfo_info[] =
{
  {
    0,
    "heap use       (MB):"
  },
  {
    1,
    "mmap use       (MB):"
  },
  {
    1,
    "sbrk use       (MB):"
  },
  {
    1,
    "malloc in use  (MB):"
  },
  {
    1,
    "malloc in free (MB):"
  }
};
# define MALLINFO_STATS sizeof mallinfo_info / sizeof mallinfo_info[0]
#endif


int ut_sysMemMallinfoNumStats(int verbosity)
{
  int numMeasurements, i;
  numMeasurements = 0;
#ifdef FLASH_SUPPORT_MALLINFO
  for (i=0; i<MALLINFO_STATS; ++i) {
    if (mallinfo_info[i].verbosity <= verbosity) {
      ++numMeasurements;
    }
  }
#endif
  return numMeasurements;
}

void ut_sysMemMallinfo(meminfo_t *meminfo, int meminfoSize, int verbosity)
{
#ifdef FLASH_SUPPORT_MALLINFO
  struct mallinfo m;
  const double convertToMB = 1.0 / (1024.0 * 1024.0);
  double memMeasurements[MALLINFO_STATS];
  int numMeasurements, i, j;

  m = mallinfo();
  memMeasurements[0] = (m.hblkhd + m.uordblks) * convertToMB; /* heap use */
  memMeasurements[1] = m.hblkhd * convertToMB; /* mmap use */
  memMeasurements[2] = m.arena * convertToMB; /* sbrk use */
  memMeasurements[3] = m.uordblks * convertToMB; /* malloc in use */
  memMeasurements[4] = m.fordblks * convertToMB; /* malloc in free */

  numMeasurements = 0;
  for (i=0; i<MALLINFO_STATS; ++i) {
    if (mallinfo_info[i].verbosity <= verbosity) {
      ++numMeasurements;
      if (numMeasurements <= meminfoSize) {
	j = numMeasurements - 1;
	meminfo[j].measurement = memMeasurements[i];
	meminfo[j].description = mallinfo_info[i].description;
      }
    }
  }
#endif
}
