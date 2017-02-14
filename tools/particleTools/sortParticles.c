#include <stdio.h>
#include <string.h>
#include "hdf5.h"
#include "mpi.h"
#include "flash_ptio.h"

#define NONEXISTENT -1
#define DOUBLE
#ifdef DOUBLE
typedef double MYREAL;
#define MPTYPE MPI_DOUBLE_PRECISION
#else
typedef float MYREAL;
#define MPTYPE MPI_REAL
#endif
#define HUGE 99999999
#define PADDING 20


int localNumParticles,maxLocalNum,numToSend,numRecv,sendCount,recvCount;
int tagUpperBound,tagLowerBound,tagIndex, numPartProps, MyPE;
int beginStep, endStep, cont;
int *tagList;
MYREAL *particlesSorted, *particlesToSend, *particlesRecd; 


void TagListAndBounds()
{
  int prevTag, i,j;
  cont=1;
  tagLowerBound=HUGE;
  tagUpperBound=0;
  j=0;
  prevTag = (int)(particlesRecd[j+tagIndex]);
  if(prevTag<tagLowerBound)
    tagLowerBound=prevTag;
  if(prevTag>tagUpperBound)
    tagUpperBound=prevTag;
  tagList[0]=prevTag;
  
  for(i=1;i<localNumParticles;i++)
    {
      j+=numPartProps;
      tagList[i]=(int)(particlesRecd[j+tagIndex]);
      if(tagList[i]<tagLowerBound)
	tagLowerBound=tagList[i];
      if(tagList[i]>tagUpperBound)
	tagUpperBound=tagList[i];       
      if(cont==1) 
	cont=tagList[i]-prevTag;
      prevTag=tagList[i];
    }
}



int findParticleIndex(int tag)
{
  int i, notFound, ind;
  if((tag<tagLowerBound) || (tag>tagUpperBound)){
    return NONEXISTENT;
  }
  else
    if(cont==1)
      ind=tag-tagLowerBound+1;
    else
      {
	ind=NONEXISTENT;
	notFound=0;
	i=0;
	while(notFound==0)
	  {
	    if(tagList[i]==tag)
	      {
		notFound=1;
		ind=i;
	      }
	    else
	      i++;
	  }  /*while*/
      } /*if*/
  return(ind-1);
}

void sortParticles(int numUnsorted)
{    
  int i,j,k,n;
  numToSend=0;
  k=0;
  for(i=0;i<numUnsorted;i++)
    {
      n = (int)(particlesRecd[k+tagIndex]);
      j=findParticleIndex(n);
      if(j != NONEXISTENT)
	{
	  j*=numPartProps;
	  for(n=0;n<numPartProps;n++)
	    particlesSorted[j+n]=particlesRecd[k+n];
	}
      else
	{
	  for(n=0; n<numPartProps; n++)
	    particlesToSend[numToSend+n]=particlesRecd[k+n];
          numToSend += numPartProps;
	}
      k += numPartProps;
    }
  if(numToSend==0)
    particlesToSend[0]=0.0;
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
  int currentFile;
  char filename[260];
  char *fileTag;

  ierr = MPI_Init(&argc,&argv);

  ierr = MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD,&MyPE);
  
  templateFile = atoi(argv[2]);
  minFile = atoi(argv[3]);
  maxFile = atoi(argv[4]);
  fileTag = argv[1];
  
  currentFile = templateFile;
  

  
  leftNegh=MyPE-1;
  rightNegh=MyPE+1;
  
  if(leftNegh<0)
    leftNegh=numProcs-1;
  if(rightNegh>=numProcs)
    rightNegh=0;

  commtag = 10;
  
  sprintf(filename, "%s%04d", fileTag, currentFile);
  flashFileOpen(filename, &file_id);  
  
  ierr = readFlashLocalNP(file_id, &localNumParticles, &maxLocalNum, 
                          &offsetStart, &offsetEnd, &tagIndex, &numPartProps, 
                          numProcs, MyPE, MPI_COMM_WORLD);
  
  if(ierr != NORMAL_STATUS){
    fprintf(stderr, "Error reading in particle metadata!\nExiting with"); 
    fprintf(stderr, " an error code of %d.\n", ierr);
    return ierr;
  }
  MPI_Allreduce(&localNumParticles, &maxLocalNum, 1, MPI_INTEGER, MPI_MAX,
                MPI_COMM_WORLD);
  
  /*getParticlesMetaData(&localNumParticles, &maxLocalNum, &tagIndex, &numPartProps);*/
  
  arraySize=1;
  if(localNumParticles>0) 
    arraySize=numPartProps*(localNumParticles+PADDING);
  
  particlesSorted = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
  
  particlesRecd = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
  particlesToSend = (MYREAL *)malloc(sizeof(MYREAL)*arraySize);
  tagList = (int *)malloc(sizeof(int) * arraySize);
  
  ierr = readFlashParticles(file_id, particlesRecd, localNumParticles, offsetStart, numPartProps);
  if(ierr != NORMAL_STATUS){
    fprintf(stderr, "Error reading in particle data!\nExiting with"); 
    fprintf(stderr, " an error code of %d.\n", ierr);
    return ierr;
  }
  TagListAndBounds();
 
  flashFileClose(&file_id);


  /*template generated.  Start sorting subsequent files*/
  currentFile = minFile;
  
  /* loop over files: first file just set our ordering */
  while(currentFile<=maxFile){
   
    for (i = 0; i < arraySize; ++i){
	particlesSorted[i] = NONEXISTENT;
    }
    
    sprintf(filename, "%s%04d", fileTag, currentFile++);
    if(MASTERPE == MyPE) 
	printf("%s\n", filename);
    flashFileOpen(filename, &file_id);        
    
    ierr = readFlashParticles(file_id, particlesRecd, localNumParticles, offsetStart, numPartProps);
    if(ierr != NORMAL_STATUS){
      fprintf(stderr, "Error reading in particle data!\nExiting with"); 
      fprintf(stderr, " an error code of %d.\n", ierr);
      return ierr;
    }
    
    sortParticles(localNumParticles);
    ierr = MPI_Allreduce(&numToSend,&maxToSend, 1, MPI_INTEGER, MPI_MAX, 
                         MPI_COMM_WORLD);
    
    
    
    count=0;
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
     
     
     writeFlashParticles(file_id, particlesSorted, localNumParticles, 
                         numPartProps, offsetStart);
     
     flashFileClose(&file_id);
     
  }/*for loop*/
  
  free(particlesSorted);
  free(particlesToSend);
  free(particlesRecd);
  free(tagList);
 
  if(MASTERPE == MyPE)
    printf("done\n");
  MPI_Finalize();
  return NORMAL_STATUS;
}







