#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include "mangle_names.h"
#include "ut_sysMem.h"

/* return the resident memory of the process in MB */
int FTOC(ut_sys_mem_usage)(double *vsizeMB, double *residentMemoryMB, int *memSampler)
{


#ifdef Darwin

  /* We can work out how to find memory usage on Macs later */
  *vsizeMB = (double) 0.0;
  *residentMemoryMB = (double) 0.0;
  *memSampler = UT_SYSMEM_NONE;

#else

#ifndef IBM /* only 2 cases: linux, and ibm.  linux does not implement getrusage, but has 
               the /proc filesystem, and ibm has a working getrusage */
  FILE *selfFile;
  
  /* data I want: virtual mem use, resident mem; 
     and stack size, heap size, static data size 
     so I know what's causing the virtual size 
     and resident use.

     For now, just get resident use.
  */

  unsigned long int vsize;
  long int rss;

  int ret;
  char s[100];
  long pagesize = sysconf(_SC_PAGESIZE);

  selfFile = fopen("/proc/self/stat", "r");
  if (!selfFile) {
    perror("io_memory_usage could not open /proc/self/stat file (io_memory_usage only works on linux, right now)");
    *residentMemoryMB = -1.0;
    *vsizeMB = -1.0;
  } else {

    ret = fscanf(selfFile,"%*d %*s %*c %*d %*d %*d %*d %*d %*lu %*lu %*lu %*lu %*lu %*lu %*lu %*ld %*ld %*ld %*ld %*ld %*ld %*lu %lu %ld %*lu ", &vsize, &rss);
    if ((ret == EOF) || (ret < 2)) {
      fprintf(stderr, "[io_memory_usage] fscanf(/proc/self/stat) returned %d!\n", ret);
    }
    
    ret = fclose(selfFile);
    if (ret != 0) {
      fprintf(stderr, "[io_memory_usage] fclose(/proc/self/stat) returned %d!\n", ret);
    }

    rss = pagesize * rss;
    /* divide by 1024 * 1024 to convert to MB */
    *residentMemoryMB = ( (double) rss )*9.765625e-4*9.765625e-4;
    *vsizeMB = ( (double) vsize )*9.765625e-4*9.765625e-4;
    *memSampler = UT_SYSMEM_PROC;
  }
#else  /* on an ibm, getrusage is supposed to work, but not sure how to get vsize */

  struct rusage r_usage;

  getrusage ( RUSAGE_SELF, &r_usage);

  *residentMemoryMB = ( (double) r_usage.ru_maxrss ) * 9.765625e-4;
  *vsizeMB = (double) 0.0;
  *memSampler = UT_SYSMEM_RUSAGE;

#endif  

#endif

  return 0;
}
