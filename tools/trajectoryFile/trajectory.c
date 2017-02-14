/* Main function for the particle trajectory file re-formater.*/
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "trajectory.h"
#include "getOptions.h"
#include "hdf5.h"
#include "mpi.h"
#include "sortParticles.h"
#include "flash_ptio.h"

extern int localNumParticles,maxLocalNum,numToSend,numRecv,sendCount,recvCount;
extern int tagUpperBound,tagLowerBound,tagIndex, blkIndex, numPartProps, MyPE;
extern int beginStep, endStep, cont;
extern int *tagList;
extern MYREAL *particlesSorted, *particlesToSend, *particlesRecd; 


int getPos(int i, int j, int k, int xdim, int ydim){

  return i + j*xdim + k*xdim*ydim;
}


int main(int argc, char **argv)
{   
  int  commtag;
  int  leftNegh,rightNegh,sendCount,recvCount;
  MPI_Status status;
  int  maxToSend;
  int  notDone, numProcs; 
  int  arraySize, ierr, count;
  int  file_id;
  int  offsetStart, offsetEnd;
  int i,j;
  int minFile, maxFile, templateFile;
  int currentFile, fileNumber;
  char filename[260];
  char *fileTag;

  hid_t inFile, outFile;
  int MyPE;

  arguments_t arguments;

  char fileNumString[10];
  char padFormat[5] = "\0\0\0\0\0";
  
  int globalNumParticles;
  
  int numFiles, offset;
  double *timesteps;
  double *outputParticles;
  int nextFlushFile, nextFlushStep, prevFlushStep, flushCount;
 
  /*for file splitting*/
  int useSplitFile = 0;
  char splitNumString[10];
  char splitNumFormat[10];
  int mySplitNum = 0;
  int splitOffset; //our offset into a split file
  int splitNumFileParticles; //the number of particles in a given split file.
  int splitPE, splitNumProcs;
  double walltime, walltime_prev, starttime;
  MPI_Comm io_comm;
  int rightSplitOffset;
  MPI_Status stat;

  /*allow us to handle an ever-shrinking number of particles*/
  int currentGlobalNumPart, currentLocalNumParticles, currentOffset;

  int step; /*keep track of how many steps we have gone through for flush purposes */

  ierr = MPI_Init(&argc,&argv);
  
  ierr = MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&MyPE);
  walltime = MPI_Wtime();
  walltime_prev = walltime;
  starttime = walltime;
  
  printf("Reading options.\n");
  /*read in and process options*/
  getOptions(&argc, argv, &arguments);


  assert(arguments.start_num <=arguments.end_num);
  assert(arguments.padDigits > 0);
  assert(arguments.splitNumber > 0);
  assert(arguments.timestepsPerFlush >= 0);
  
  
  numFiles = arguments.end_num - arguments.start_num + 1;
  
  if(MASTERPE == MyPE)
    printf("numFiles = %d\n", numFiles);

  /*file splitting turned on*/
  io_comm = MPI_COMM_WORLD;
  if(arguments.splitNumber > 1){
    printf("file splitting enabled for %d files\n.",arguments.splitNumber);
    useSplitFile = 1;
    strncpy(splitNumFormat, "_s%05d", 10);
    
    mySplitNum = MyPE * arguments.splitNumber/numProcs; //+
      //(MyPE > (numProcs/arguments.splitNumber) ? 1:0);
    //printf("MyPE: %d : Using spliting.  My number: %d\n",MyPE, mySplitNum);
    sprintf(splitNumString, splitNumFormat, mySplitNum);
    /*use same ordering as the default communicator*/
    MPI_Comm_split(MPI_COMM_WORLD, mySplitNum, 0, &io_comm);
    MPI_Comm_size(io_comm,&splitNumProcs);
    MPI_Comm_rank(io_comm,&splitPE);  
  }
  else{
    splitNumProcs = numProcs;
    splitPE = MyPE;
  }

  
  /*open output file*/
  sprintf(padFormat, "%%0%dd", arguments.padDigits);
  strncat(strncpy(filename, arguments.basename,256),"traj",260);
  if(arguments.splitNumber > 1)
    strncat(filename, splitNumString, 260);
  flashFileOpen(filename, &outFile, 1, io_comm);
  
  

  


  /*neighbor determination*/
  leftNegh=MyPE-1;
  rightNegh=MyPE+1;
  
  if(leftNegh<0)
    leftNegh=numProcs-1;
  if(rightNegh>=numProcs)
    rightNegh=0;

  commtag = 10;
  
  /*set the initial file to flush out the buffers on*/
  walltime = MPI_Wtime();
  printf("beginning loop over files, took %f seconds.\n", walltime - walltime_prev);
  walltime_prev = walltime;

  step = 0;
  flushCount = 0;
  /*loop over files, template off of the first file*/
  for(fileNumber = arguments.start_num; 
      fileNumber <= arguments.end_num; 
      fileNumber++){
    

    if(MASTERPE == MyPE)
      printf("Begin processing file number %d.\n", fileNumber);

    /*open input file*/
    sprintf(fileNumString, padFormat, fileNumber);
    strncat(strncpy(filename, arguments.basename,256),fileNumString,260);
    flashFileOpenReadOnly(filename, &inFile);
    
    if(fileNumber == arguments.start_num){ /*only for first file*/

      copyUnknownNames(&inFile, &outFile, splitPE, &tagIndex, & blkIndex);
      

      MPI_Bcast(&tagIndex, 1, MPI_INT, MASTERPE, MPI_COMM_WORLD);
      MPI_Bcast(&blkIndex, 1, MPI_INT, MASTERPE, MPI_COMM_WORLD);
      assert(tagIndex >= 0);
      assert(blkIndex >= 0);
      
      /*this number does change through a run.  Can decrease.*/
      getNumberOfParticles(&inFile, &globalNumParticles, &numPartProps);
      if(MASTERPE == MyPE){
	printf("numParticles %d\n", globalNumParticles);
	printf("numPartProps %d\n", numPartProps);
      }

      /*local number of particles and offsets*/
      currentLocalNumParticles = localNumParticles = (globalNumParticles / numProcs) + (((globalNumParticles % numProcs) > MyPE) ? 1:0);
      //printf("lnp: %d\n", localNumParticles);
      if(arguments.splitNumber > 1){
	MPI_Allreduce(&localNumParticles, &splitNumFileParticles, 1,
		      MPI_INT, MPI_SUM, io_comm);
	//printf("num parts in pe %d file: %d\n", MyPE, splitNumFileParticles);
      }
      else 
	splitNumFileParticles = globalNumParticles;
   
      /*we can calculate the number to the left from this decomposition*/
      currentOffset = offset = MyPE*(globalNumParticles/numProcs) + 
        (((globalNumParticles % numProcs) > MyPE)? 
         MyPE:(globalNumParticles % numProcs));
      
      /*find the split offset of my splitPE*/
      /*this is the offset into the split file.*/
      
      /*pass to the right neighbor on your communicator.*/

      if(splitPE == 0)
	splitOffset = 0;
      
      if(splitPE > 0){
	/*get from left*/
	MPI_Recv(&splitOffset,1,MPI_INT,splitPE-1,1,
		 io_comm,&stat);
      }
      if(splitPE < splitNumProcs - 1){
	/*send to right*/
	rightSplitOffset = splitOffset + localNumParticles;
	MPI_Send(&rightSplitOffset,1,MPI_INT,splitPE+1,1,
		 io_comm);
      }
	
      //printf("%d.%d: splitOffset = %d\n",MyPE,splitPE, splitOffset);
      
      
      //This doesn't work --PMR
      //splitOffset = splitPE*(splitNumFileParticles/splitNumProcs) + 
      //(((splitNumFileParticles % splitNumProcs) > splitPE)? 
      // splitPE:(splitNumFileParticles % splitNumProcs));
      
      //printf("%d: offset: %d\n",MyPE, offset);
      //printf("%d: splitOffset: %d\n",MyPE, splitOffset);
      

      /*initialize the trajectory dataset in the output file*/
      trajectoryInit(outFile, numPartProps, numFiles, splitNumFileParticles); 
      
      /*allocate bufferes for file iteration*/
      arraySize=1;
      if(localNumParticles>0) 
	arraySize=numPartProps*(localNumParticles+PADDING);
      
      particlesSorted = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
      particlesRecd = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
      particlesToSend = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
      tagList = (int *)malloc(sizeof(int) * arraySize);
      
      assert(particlesSorted != NULL);
      assert(particlesRecd != NULL);
      assert(particlesToSend != NULL);
      assert(tagList != NULL);
      
      /*limit the number files we hold at once.  We flush every n file reads*/
      if(0 >= arguments.timestepsPerFlush || arguments.timestepsPerFlush > numFiles)
	arguments.timestepsPerFlush = numFiles;
      
      // printf("TimestepsPerFlush, %d\n",arguments.timestepsPerFlush);
      nextFlushStep = arguments.timestepsPerFlush;
      prevFlushStep = 0;

      outputParticles = (double *)malloc(sizeof(double)*numPartProps*arguments.timestepsPerFlush*localNumParticles);
      timesteps = (double *)malloc(sizeof(double)*arguments.timestepsPerFlush);
      assert(outputParticles != NULL);
      assert(timesteps != NULL);
      /*generate a list of tags -- These are from the first file*/
      
      readFlashParticles(inFile, particlesRecd, localNumParticles, offset, 
			 numPartProps);
      TagListAndBounds();
      
      /*dump out particles from this file.  By definition, these are ordered.*/
      /*these go to the marshalling array*/
		      
      for( i = 0; i < localNumParticles; i++){
      	memcpy(&outputParticles[numPartProps*(flushCount) + i*numPartProps*arguments.timestepsPerFlush],
             &particlesRecd[i*numPartProps],
             sizeof(double)*numPartProps);
	
      }
      

      /*copy our to sorted for step n+1*/
      memcpy(particlesSorted, particlesRecd, sizeof(double)*numPartProps*localNumParticles);

      
      getTimeStep(&inFile, &timesteps[fileNumber - arguments.start_num]);
    
      flashFileClose(&inFile);
   
    }/*end if(fileNumber == arguments.start_num)*/
   

    else{ /* we are no longer on the first iteration */

    /*set blk to -1.  This will allow us to easily see what has left the 
      domain, while preserving tag.*/
    for (i = 0; i < localNumParticles; ++i){
      particlesSorted[blkIndex + i*numPartProps] = NONEXISTENT;
    } 

    /*has our globalNumberOfParticles changed, if so recompute.*/
    getNumberOfParticles(&inFile, &currentGlobalNumPart, &numPartProps);
    
    /*copy n-1 into rec array */
    /*we've changed, need to recompute*/
   
    if(currentGlobalNumPart != globalNumParticles){
      globalNumParticles = currentGlobalNumPart;
      
      currentLocalNumParticles = (globalNumParticles / numProcs) + 
	(((globalNumParticles % numProcs) > MyPE) ? 1:0);
      
      currentOffset = MyPE*(globalNumParticles/numProcs) + 
        (((globalNumParticles % numProcs) > MyPE)? 
         MyPE:(globalNumParticles % numProcs));
    }
      
    /*gather and sort particles*/
    ierr = readFlashParticles(inFile, particlesRecd, currentLocalNumParticles,
			      currentOffset, numPartProps);
    if(ierr != NORMAL_STATUS){
      fprintf(stderr, "Error reading in particle data!\nExiting with"); 
      fprintf(stderr, " an error code of %d.\n", ierr);
      return ierr;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("prior to entry MyPe = %d\n", MyPE);
    sortParticles(currentLocalNumParticles);
    
    ierr = MPI_Allreduce(&numToSend,&maxToSend, 1, MPI_INTEGER, MPI_MAX, 
                         MPI_COMM_WORLD);
    count=0;
    //printf("%d: all reduce complete\n", MyPE);
    
    
    /*sieve sorting*/
    while (maxToSend>0)
      {
        recvCount = (maxToSend)*numPartProps;
	if(numToSend>0)
	  sendCount = numToSend;
	else
	  sendCount = 1;
	
	ierr = MPI_Sendrecv(particlesToSend,sendCount,MPTYPE,rightNegh,commtag,
			    particlesRecd,recvCount,MPTYPE,leftNegh,commtag,
			    MPI_COMM_WORLD, &status);
	ierr =  MPI_Get_count(&status,MPTYPE,&recvCount);
	numRecv=recvCount/numPartProps;
	
	if(numRecv>=0)
	  sortParticles(numRecv);
	ierr = MPI_Allreduce(&numToSend,&maxToSend, 1, MPI_INTEGER, MPI_MAX, 
			     MPI_COMM_WORLD);
	
	count++;
	if(count>numProcs)
	  {
	    for(j=tagIndex; j<numToSend; j+=numPartProps)
	      printf("attempting to send part %f on proc %d\n", particlesToSend[j], MyPE);
	    printf("%d to send, should be 0! pe: %d\n", numToSend/numPartProps, MyPE);
	    fflush(stdout);
	    MPI_Barrier(MPI_COMM_WORLD);
	    MPI_Abort(MPI_COMM_WORLD, -1);
	  }	
      }
    
    for( i = 0; i < localNumParticles; i++){
      memcpy(&outputParticles[numPartProps*(flushCount) + i*numPartProps*(nextFlushStep - prevFlushStep)],
             &particlesSorted[i*numPartProps],
             sizeof(double)*numPartProps);
      
    }


    walltime = MPI_Wtime();
    printf("File Number %d sorted, took %f seconds.\n",fileNumber ,walltime-walltime_prev);
    walltime_prev = walltime;


    /*particlesSorted implicitly will have step n-1 on the next iteration*/
    
    
    /*extract timestep*/
    getTimeStep(&inFile, &(timesteps[flushCount]));
    
    flashFileClose(&inFile);

    }/*end else for first-step if*/
 
    flushCount++;
    step++;
    
    /*DEV: this is for debugging*/
    if (step >= nextFlushStep){
      //if(MASTERPE == MyPE)
	//printf ("nextFlushStep = %d, prevFlushStep = %d, numToWrite: %d\n", nextFlushStep, prevFlushStep, nextFlushStep - prevFlushStep);
      
      for (i = 0; i < localNumParticles; ++i){
	for(j = 0; j < nextFlushStep - prevFlushStep; ++j){
	  assert(outputParticles[tagIndex + j*numPartProps + i*numPartProps*arguments.timestepsPerFlush] > 0);
	}
      } 
      assert(localNumParticles > 0);
      if(splitPE == 0)
	assert(splitOffset == 0);
      
      //printf ("writing from PE: %d, splitOffset = %d, localNumParticles= %d\n", MyPE, splitOffset, localNumParticles);
      writeTrajectoryData(&outFile, outputParticles, splitOffset, 
			  numPartProps, prevFlushStep, 
			  nextFlushStep - prevFlushStep, arguments.timestepsPerFlush, 
			  splitNumFileParticles, 
			  localNumParticles);
      //memset(outputParticles, 0, sizeof(double)*arguments.timestepsPerFlush*localNumParticles*numPartProps);
      writeTimestepData(outFile, nextFlushStep - prevFlushStep, prevFlushStep, timesteps);
      //printf("Step: %d\n", step);
      prevFlushStep = nextFlushStep;
      nextFlushStep += arguments.timestepsPerFlush;
      if(nextFlushStep > numFiles){
	//prevFlushStep++;
	nextFlushStep = numFiles;
      }
      flushCount = 0;
      // printf ("Post write: nextFlushStep = %d, prevFlushStep = %d\n", nextFlushStep, prevFlushStep);
      //break;
    }
  } /* end loop over files */  
  
 
  //writeTimestepData(&outFile, &numFiles, timesteps);

  

  free(outputParticles);
  free(tagList);
  free(particlesToSend);
  free(particlesSorted);
  free(particlesRecd);

  MPI_Barrier(io_comm);
  //printf("attempring closes\n");
  flashFileClose(&outFile);
  //printf("close complete\n");

  if(MASTERPE == MyPE){
    walltime = MPI_Wtime();
    printf("done, runtime was %f seconds\n", walltime - starttime);
  }
  
  MPI_Finalize();
  
  return NORMAL_STATUS;
}

