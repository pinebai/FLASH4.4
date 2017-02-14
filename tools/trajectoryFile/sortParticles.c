#include <stdio.h>
#include <string.h>
#include "hdf5.h"
#include "mpi.h"
#include "flash_ptio.h"
#include "sortParticles.h"




int localNumParticles,maxLocalNum,numToSend,numRecv,sendCount,recvCount;
int tagUpperBound,tagLowerBound,tagIndex, blkIndex, numPartProps, MyPE;
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
  int i, notFound, ind,j;
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
	for(j=0;j<localNumParticles;++j){
	  if(tagList[j] == tag)
	    return i;
	}
	return NONEXISTENT;
	//while(notFound==0)
	//{
	//  if(tagList[i]==tag)
	//    {
	//notFound=1;
	//ind=i;
	//return ind;
	//    }
	//  else
	//    i++;
	// }  /*while*/
      } /*if*/
  return(ind-1);
}
#include<stdlib.h>

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











