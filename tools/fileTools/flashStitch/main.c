#include "options.h"
#include "flashStitcher.h"

int parse_cmdline(options_t * opts, int argc, char **argv);

int main(int argc, char **argv){
  int retValue = 0;

  options_t opts;
  
  retValue = parse_cmdline(&opts, argc, argv);
  if (retValue < 0){
    return retValue;
  }

#ifdef DEBUG
  printf("%s %d %d\ndone\n", opts.filename, opts.splitNum, opts.zeroPad);
#endif

   
  retValue = 
    flashStitcher(opts.filename, opts.splitNum,  opts.zeroPad, &opts);
  
  if(retValue != 0){
    printf("Abnormal exit condition!\n");
  }
  
  return retValue;
}

