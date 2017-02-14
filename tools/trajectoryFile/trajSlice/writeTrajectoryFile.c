#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "hdf5.h"
#include "flash_ptio.h"

//extern int MyPE;

int fileOpenRead(char *filename, hid_t *fid){

  *fid = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
  assert(*fid > 0);

  return NORMAL_STATUS;
}

int fileOpenWrite(char *filename, hid_t *fid){

  *fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  assert(*fid > 0);
  return NORMAL_STATUS;
}


int fileClose(hid_t *fid){

  return H5Fclose(*fid);
  
}



int copyUnknownNames(hid_t fidIn, hid_t fidOut, int* tagIndex, int* blkIndex,
		     int useThreshold, char *thresholdName, 
		     int *thresholdIndex){

  hsize_t dimens[2], maxdimens[2], stride[2], count[2], start[2];
  hid_t dataspace, dataset, memspace, string_type, plist_id;

  int rank, i;
  int ierr;

  char *partNames;

  *tagIndex = -1;
  *blkIndex = -1;
  *thresholdIndex = -1;

  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, 24);

  dataset = H5Dopen(fidIn, "particle names");
  
  if(dataset < 0){
    printf("dataset open failed!\n");
    return DATASET_OPEN_FAIL;
  }
  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens, maxdimens);
  
  rank = 2;
  
  start[0] = (hsize_t)0;
  start[1] = (hsize_t)0;
  
  stride[0] = (hsize_t)1;
  stride[1] = (hsize_t)1;
  
  count[0] = dimens[0];
  count[1] = 1;
  
  partNames = (char *)malloc(sizeof(char)*(int)dimens[0]*24);
  
  
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);
  memspace = H5Screate_simple(rank, dimens, NULL);
  
  
  H5Dread(dataset, string_type, memspace,dataspace, H5P_DEFAULT, partNames);
  
  
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  /*get tag index*/
  for(i = 0; i < count[0]; i++){
    if(!(strncmp(&partNames[i*24], "tag ", 4))){
      *tagIndex = i;
    }
    else if(!(strncmp(&partNames[i*24], "blk ", 4))){
      *blkIndex = i;
    }
    else if(useThreshold && (!(strncmp(&partNames[i*24], thresholdName, 4)))){
      *thresholdIndex = i;
    }
  }
  
  
  /*write back out to trajectory file*/
  
  dataspace = H5Screate_simple(rank, dimens, NULL);
  if(dataspace < 0)
    printf("dataspace error\n");
  
  dataset = H5Dcreate(fidOut, "particle names", string_type, dataspace, H5P_DEFAULT);
  if(dataset < 0){
    printf("Could not open dataset for writing!\n");
  }
  

  i = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, partNames);
  if (i < 0)
    printf("write failed!");
  
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(string_type);
  
  free(partNames);
  
  return NORMAL_STATUS;
}


int copyTimesteps(hid_t fidIn, hid_t fidOut){

  hid_t memspace, dataspace, dataset;
  hsize_t dimens, start, count, stride;
  
  double *timesteps;

  dataset = H5Dopen(fidIn, "timesteps");
  dataspace = H5Dget_space(dataset);

  H5Sget_simple_extent_dims(dataspace, &dimens, NULL);

  timesteps = (double *)malloc(sizeof(double)*dimens);
  
  memspace = H5Screate_simple(1, &dimens, NULL);

  H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
	     H5P_DEFAULT, timesteps);
  
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memspace);

  /*write stage*/

  dataspace = H5Screate_simple(1, &dimens, NULL);
  dataset = H5Dcreate(fidOut, "timesteps", H5T_NATIVE_DOUBLE,
		      dataspace, H5P_DEFAULT);
 
  
  memspace = H5Screate_simple(1, &dimens, NULL);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
	     H5P_DEFAULT, timesteps);
  
  
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

 

  free(timesteps);

  return NORMAL_STATUS;
}


int copyParticleSelection(hid_t fidIn, hid_t fidOut, 
			  int partStart, int partStride, int partCount){

  hid_t dataspace, memspace, dataset;
  hsize_t dimens[3], start[3], stride[3], count[3];

  double *trajectories;
  /*read in a set of particles that we are interested in*/

  
  dataset = H5Dopen(fidIn, "particle trajectories");
  dataspace = H5Dget_space(dataset);

  H5Sget_simple_extent_dims(dataspace, dimens, NULL);

  start[0] = (hsize_t)partStart;
  start[1] = 0;
  start[2] = 0;
  
  stride[0] = (hsize_t)partStride;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = (hsize_t)partCount;
  count[1] = dimens[1];
  count[2] = dimens[2];

  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);
  
  trajectories = (double *)malloc(sizeof(double)*count[0]*count[1]*count[2]);

  memspace = H5Screate_simple(3, count, NULL);

  H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, 
	  trajectories);

  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Dclose(dataset);
 

  /*write stage*/

  dataspace = H5Screate_simple(3, count, NULL);
  dataset = H5Dcreate(fidOut, "particle trajectories", H5T_NATIVE_DOUBLE,
		      dataspace, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  	   trajectories);
  
  H5Dclose(dataset);
  H5Sclose(dataspace);
 

  free(trajectories);

  return NORMAL_STATUS;
}


/* Copies particles based on a threshold value.  Right now this is strictly
   less than.  Can add comparitors later.
   Arguments:
   fidIn: file to read
   fidOut: file to write to
   maxCount: an integer number of particles to copy
   maxThreshold: an integer number of particles to apply the threshold 
                 condition to.  If this is smaller than maxCount, the
		 function just copies the next maxCount - maxThreshold particles
   thresholdIndex: Index of the variable to apply the threshold to.
   thresholdValue: if the trajectory ever has a value > than this it is coppied
                   by the threshold condition, up to maxCount or maxThreshold, 
		   whichever comes first.
   Returns: 0 on a successful execution.
*/

int copyParticleThreshold(hid_t fidIn, hid_t fidOut, 
			  int maxCount, int maxThreshold, int thresholdIndex, 
			  double thresholdValue){

  hid_t dataspaceRead, memspaceRead, datasetRead;
  hid_t dataspaceWrite, memspaceWrite, datasetWrite;
  hsize_t dimens[3], start[3], stride[3], count[3];

  int i, partsWritten, currentPart, maxPartsInFile;
  
  double *trajectory;
  /*read in a set of particles that we are interested in*/

  
  datasetRead = H5Dopen(fidIn, "particle trajectories");
  dataspaceRead = H5Dget_space(datasetRead);
  
  
  H5Sget_simple_extent_dims(dataspaceRead, dimens, NULL);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  
  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = maxCount;
  count[1] = dimens[1];
  count[2] = dimens[2];

  maxPartsInFile = dimens[0];

  trajectory = (double *)malloc(sizeof(double)*count[1]*count[2]);

  memspaceRead = H5Screate_simple(2, &dimens[1], NULL);
 
  /*write stage*/
  
  dataspaceWrite = H5Screate_simple(3, count, NULL);
  datasetWrite = H5Dcreate(fidOut, "particle trajectories", H5T_NATIVE_DOUBLE,
		      dataspaceWrite, H5P_DEFAULT);
  

  //reopen for write process

  currentPart = partsWritten = 0;
  count[0] = 1;
  while (partsWritten < maxThreshold && currentPart < maxPartsInFile && 
	 partsWritten < maxCount){

    //read
    start[0] = currentPart;
    start[1] = 0;
    start[2] = 0;

    H5Sselect_hyperslab(dataspaceRead, H5S_SELECT_SET, 
			    start, stride, count, NULL);

    H5Dread(datasetRead, H5T_NATIVE_DOUBLE, memspaceRead, dataspaceRead, 
	    H5P_DEFAULT, trajectory);
    //compare
    
    for(i = 0; i < count[1]; i++){
    // printf("grabbing element %d\n", (thresholdIndex+ i*(int)count[2]));
    if (trajectory[thresholdIndex + i*(int)count[2]] >
      thresholdValue){
	
      //write
	start[0] = partsWritten;
	
	H5Sselect_hyperslab(dataspaceWrite, H5S_SELECT_SET, 
			    start, stride, count, NULL);
  
	H5Dwrite(datasetWrite, H5T_NATIVE_DOUBLE, memspaceRead, dataspaceWrite, 
		 H5P_DEFAULT, trajectory);

	
    partsWritten++;
    break;
    }
  }
    currentPart++;
  }

  /* this is so we can populate the rest of the file.*/
  if (maxCount > maxThreshold){

    printf("Threshold particles filled.\n");
    while(partsWritten < maxCount && currentPart < maxPartsInFile){

      /*just copy the next trajectory until full, or we run out of 
        trajectories*/
      //read
      start[0] = currentPart;
      start[1] = 0;
      start[2] = 0;
      
      H5Sselect_hyperslab(dataspaceRead, H5S_SELECT_SET, 
			  start, stride, count, NULL);
      
      H5Dread(datasetRead, H5T_NATIVE_DOUBLE, memspaceRead, dataspaceRead, 
	      H5P_DEFAULT, trajectory);
      
	
      //write
      start[0] = partsWritten;
      
      H5Sselect_hyperslab(dataspaceWrite, H5S_SELECT_SET, 
			  start, stride, count, NULL);
      
      H5Dwrite(datasetWrite, H5T_NATIVE_DOUBLE, memspaceRead, dataspaceWrite, 
	       H5P_DEFAULT, trajectory);

	
      partsWritten++;
      currentPart++;
      
    }

  }

   
  H5Sclose(dataspaceRead);
  H5Sclose(memspaceRead);
  H5Dclose(datasetRead);
 
  H5Dclose(datasetWrite);
  H5Sclose(dataspaceWrite);

  free(trajectory);

  return NORMAL_STATUS;
}

int compareGreater(double a, double b){
  return (a > b)? 1 : 0;
}

int compareLess(double a, double b){
  return (a < b)? 1 : 0;
}

/* returns a count of threshold trajectories
   This allows us to size a dataset correctly for extracting all of the
   "deflagration" trajectories.
   fid : fileID to read
   thresholdIndex: an integer index for the appropriate variable to compare.
   thresholdValue: our value to compare against
   count: an integer number of trajectories that conform to the comparison.
          zeroed on each call.
   compare: a function pointer for easily changing greater/less comparisons.

   returns 0 on successful execution.
*/

int findThresholdTrajectories(hid_t fid, int thresholdIndex, 
			      double thresholdValue, 
			      int *partCount, int(*compare)(double, double)){
  
  hsize_t start[3], stride[3], dimens[3], count[3];
  hid_t memspace, dataspace, dataset;
  
  double *trajectory;
  int i, partsRead;
  
  *partCount = 0;
  
  dataset = H5Dopen(fid, "particle trajectories");
  dataspace = H5Dget_space(dataset);

  H5Sget_simple_extent_dims(dataspace, dimens, NULL);

  trajectory = (double *)malloc(sizeof(double)*dimens[1]*dimens[2]);

  memspace = H5Screate_simple(2, &dimens[1], NULL);


  start[0] = 0;
  start[1] = 0;
  start[2] = 0;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = 1;
  count[1] = dimens[1];
  count[2] = dimens[2];

  partsRead = 0;
  

  while(partsRead < dimens[0]){
    
    start[0] = partsRead;
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			start, stride, count, NULL);
    
    H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,
	    trajectory);

    for(i = 0; i < count[1]; i++){

      if (compare(trajectory[thresholdIndex + i*count[2]],
		  thresholdValue)){
	(*partCount)++;
	printf("def found: %d\n", *partCount);
	break;
      }
      
    }
    partsRead++;
    printf("part %d read\n", partsRead);
  }

  printf("Parts read: %d Parts to xfer %d\n", partsRead, *partCount);

  
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
    
  free(trajectory);
  
  return 0;
}
