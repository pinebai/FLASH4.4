#include "ut_sysMemProc.h"

static const char * proc_description[] =
{
  "/proc rss      (MB):",
  "/proc vsize    (MB):"
};

int ut_sysMemProcNumStats(int verbosity)
{
  return PROC_STATS;
}

void ut_sysMemProc(meminfo_t *meminfo, int meminfoSize, int verbosity)
{
#ifdef FLASH_SUPPORT_PROC
  const double convertToMB = 1.0 / (1024.0 * 1024.0);
  double memMeasurements[PROC_STATS];
  FILE *selfFile;
  unsigned long int vsize;
  long int rss;
  int ret, numMeasurements, i;

  long pagesize = sysconf(_SC_PAGESIZE);

  selfFile = fopen("/proc/self/stat", "r");
  if (selfFile) {

    ret = fscanf(selfFile,"%*d %*s %*c %*d %*d %*d %*d %*d %*lu %*lu %*lu %*lu %*lu %*lu %*lu %*ld %*ld %*ld %*ld %*ld %*ld %*lu %lu %ld %*lu ", &vsize, &rss);
    if ((ret == EOF) || (ret < 2)) {
      fprintf(stderr, "[ut_sysMemProc] fscanf(/proc/self/stat) returned %d!\n", ret);
    }

    ret = fclose(selfFile);
    if (ret != 0) {
      fprintf(stderr, "[ut_sysMemProc] fclose(/proc/self/stat) returned %d!\n", ret);
    }

    rss = pagesize * rss;
    memMeasurements[0] = ( (double) rss ) * convertToMB;
    memMeasurements[1] = ( (double) vsize ) * convertToMB;

    numMeasurements = meminfoSize < PROC_STATS ? meminfoSize : PROC_STATS;
    for (i=0; i<numMeasurements; ++i) {
      meminfo[i].measurement = memMeasurements[i];
      meminfo[i].description = proc_description[i];
    }
  }
#endif
}
