#include <stdlib.h>
#include <stdio.h>
#include "hdf5.h"
#include "flash_ptio.h"
#include "getOptions.h"

int main (int argc, char **argv){


  int numRead;
  hid_t fileIn, fileOut;
  int tagIndex, blkIndex, thresholdIndex;
  struct arguments_t args;
  
  
  getOptions(&argc, argv, &args);

  /*open input and output files*/
  
  fileOpenRead(args.filenameIn, &fileIn);
  fileOpenWrite(args.filenameOut, &fileOut);

  copyUnknownNames(fileIn, fileOut, &tagIndex, &blkIndex, args.useThreshold,
		   args.thresholdVarName, &thresholdIndex);
  printf("copied particle names.\n");

  copyTimesteps(fileIn, fileOut);
  printf("copied timesteps.\n");

  if (args.useThreshold){
    printf("starting threshold copy\n");

    if(args.findAllThreshold){
      findThresholdTrajectories(fileIn, thresholdIndex, 0.1, 
				&args.numToRead, compareGreater);
      args.start = args.numToRead;
    }

    printf("start == %d, read == %d\n",args.start,args.numToRead);
    //args.start is overridden to provide the maximum parts to provide a
    //threshold to.

    copyParticleThreshold(fileIn, fileOut, args.numToRead, args.start, 
			  thresholdIndex, 0.1);
 
  }
  else{
    copyParticleSelection(fileIn, fileOut, args.start, args.stride, args.numToRead);
  }
  printf("copied selection.\n");
		
  fileClose(&fileIn);
  fileClose(&fileOut);

  printf("done.\n");

}
