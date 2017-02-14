#include "ut_sysMemRusage.h"

static const char * rusage_description[] =
{
  "rusage rss     (MB):"
};

int ut_sysMemRusageNumStats(int verbosity)
{
  return RUSAGE_STATS;
}

void ut_sysMemRusage(meminfo_t *meminfo, int meminfoSize, int verbosity)
{
#ifdef FLASH_SUPPORT_RUSAGE
  const double convertToMB = 1.0 / 1024.0;
  double memMeasurements[RUSAGE_STATS];
  int numMeasurements, i;
  struct rusage r_usage;

  getrusage (RUSAGE_SELF, &r_usage);
  /* The measurement is in kilobytes */
  memMeasurements[0] = ( (double) r_usage.ru_maxrss ) * convertToMB;
  numMeasurements = meminfoSize < RUSAGE_STATS ? meminfoSize : RUSAGE_STATS;
  for (i=0; i<numMeasurements; ++i) {
    meminfo[i].measurement = memMeasurements[i];
    meminfo[i].description = rusage_description[i];
  }
#endif
}
