#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hdf5.h"
#include "mpi.h"
#include "flash_ptio.h"

extern int MyPE;

int copyUnknownNames(hid_t *fidIn, hid_t *fidOut, int myPE, int* tagIndex, int* blkIndex){

  hsize_t dimens[2], maxdimens[2], stride[2], count[2], start[2];
  hid_t dataspace, dataset, memspace, string_type, plist_id;

  int rank, i;
  int ierr;

  char *partNames;

  *tagIndex = -1;
  *blkIndex = -1;
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, 24);

  dataset = H5Dopen(*fidIn, "particle names");
  
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
  
  if(myPE == MASTERPE){
    H5Dread(dataset, string_type, memspace,dataspace, H5P_DEFAULT, partNames);
  }
  
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
  }
  
  
  /*write back out to trajectory file*/
  
  dataspace = H5Screate_simple(rank, dimens, NULL);
  if(dataspace < 0)
    printf("dataspace error\n");
  
  dataset = H5Dcreate(*fidOut, "particle names", string_type, dataspace, H5P_DEFAULT);
  if(dataset < 0){
    printf("Could not open dataset for writing!\n");
  }
  
  if(MASTERPE == myPE){
    i = H5Dwrite(dataset, string_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, partNames);
    if (i < 0)
      printf("write failed!");
  }
  
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Tclose(string_type);
  
  free(partNames);
  
  return NORMAL_STATUS;
}

/*extract the current timestep from a particle file*/
int getTimeStep(hid_t *fidIn, double *timestep){

  hid_t dataspace, memspace, dataset, string_type, real_list_type;
  hsize_t dimens, maxdimens, count, stride, step;
  int i;
  real_list_t *list;
  
  string_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(string_type, LIST_STRING_SIZE);

  dataset = H5Dopen(*fidIn, "real scalars");
  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, &dimens,&maxdimens);
  
  list = (real_list_t *)malloc(sizeof(real_list_t)*dimens);

  memspace = H5Screate_simple(1, &dimens, NULL);
  
  /*type definition for HDF5*/
  real_list_type = H5Tcreate(H5T_COMPOUND, sizeof(real_list_t));
  
  H5Tinsert(real_list_type, 
	    "name", 
	    HOFFSET(real_list_t, name), 
	    string_type);
  H5Tinsert(real_list_type, 
	    "value", 
	    HOFFSET(real_list_t, value), 
	    H5T_NATIVE_DOUBLE);
  if(MASTERPE == MyPE){
    H5Dread(dataset, real_list_type, memspace, dataspace,
		   H5P_DEFAULT, list);
  }

  for (i = 0; i < dimens; i++){
    if(!(strncmp(list[i].name, "time", 4))){
	*timestep = list[i].value;
	break;
      }
  
  }

  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memspace);
  H5Tclose(string_type);
  H5Tclose(real_list_type);
  free(list);

  return NORMAL_STATUS;
}

int writeTimestepData(hid_t fidOut, int numTimeSteps, int offset, double* timesteps){

  hid_t memspace, dataspace, dataset;
  
  hsize_t dimens, start, count, stride;
  int i;

  dimens = numTimeSteps;

  dataset = H5Dopen(fidOut, "timesteps");
  dataspace = H5Dget_space(dataset);

  start = (hsize_t)offset;
  stride = 1;
  count = (hsize_t)numTimeSteps;
  
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start, &stride, &count, NULL);
  
  start = 0;
  
  memspace = H5Screate_simple(1, &dimens, NULL);
  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, &start, &stride, &count, NULL);

  /*proc 0 really should be the only thing writing this */
  if(MASTERPE == MyPE){
    //  for(i = 0; i < numTimeSteps; i++)
    //printf("outputing timesteps: %e\n", timesteps[i]);
    
    H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
	     H5P_DEFAULT, timesteps);
  }
  
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  
  return NORMAL_STATUS;
}

int writeTrajectoryData(hid_t *fidOut, double *particleData, int offset, int numPartProps, int timestepStart, int timeSteps, int defaultStride, int globalNumParticles, int localNumParticles){
  
  hid_t memspace, dataspace, dataset;
  int totTimesteps;
  hsize_t start[3], stride[3], count[3], dimens[3];
  static int firstCall = 1;
  
  //dimens[2] = (hsize_t) numPartProps;
  //dimens[1] = (hsize_t) timeSteps;
  //dimens[0] = (hsize_t) globalNumParticles;

  dataset = H5Dopen(*fidOut, "particle trajectories");
  dataspace = H5Dget_space(dataset);

  //H5Sget_simple_extent_dims(dataspace, dimens, NULL);
  //totTimesteps =(int)dimens[1];

  /*select from dataspace*/
  start[0] = (hsize_t) offset;
  start[1] = timestepStart;
  start[2] = 0;

  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = (hsize_t) localNumParticles;
  count[1] = (hsize_t) timeSteps;
  count[2] = (hsize_t) numPartProps;

  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, stride, count, NULL);

  /*select from memory space*/
  dimens[0] = (hsize_t) localNumParticles;
  dimens[1] = (hsize_t) timeSteps;
  dimens[2] = (hsize_t) numPartProps;
  memspace = H5Screate_simple(3, dimens, NULL);

  start[0] = 0;
  start[1] = 0;
  start[2] = 0;
  
  stride[0] = 1;
  stride[1] = 1;
  stride[2] = 1;

  count[0] = (hsize_t) localNumParticles;
  count[1] = (hsize_t) timeSteps;
  count[2] = (hsize_t) numPartProps;



  H5Sselect_hyperslab(memspace, H5S_SELECT_SET, start, stride, count, NULL);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace,
	   H5P_DEFAULT, particleData);
  
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return NORMAL_STATUS;
}


int trajectoryInit(hid_t fid, int numPartProps, int timesteps, int globalNumParticles){

  hid_t dataset, dataspace;
  hsize_t dimens[3];
  int rank;
  
  rank = 3;
  dimens[0] = (hsize_t) globalNumParticles;
  dimens[1] = (hsize_t) timesteps;
  dimens[2] = (hsize_t) numPartProps;
  
  dataspace = H5Screate_simple(3, dimens, NULL);

  dataset = H5Dcreate(fid, "particle trajectories", H5T_NATIVE_DOUBLE, 
			dataspace, H5P_DEFAULT);  

  H5Dclose(dataset);
  H5Sclose(dataspace);

  /*init the timestep array, too*/
  
  dataspace = H5Screate_simple(1, &(dimens[1]), NULL);
  dataset = H5Dcreate(fid, "timesteps", H5T_NATIVE_DOUBLE, 
		      dataspace, H5P_DEFAULT);
  
  H5Dclose(dataset);
  H5Sclose(dataspace);

  return NORMAL_STATUS;
}

