#ifndef FLASH_PTIO_H
#define FLASH_PTIO_H


#define MASTERPE 0

/*error codes*/
#define NORMAL_STATUS 0
#define DATASET_OPEN_FAIL -1
#define FILE_OPEN_FAIL -2
#define HYPERSLAB_SELECT_FAIL -3
#define MEMSPACE_SELECT_FAIL -4
#define DATASPACE_SELECT_FAIL -5
#define DATA_READ_ERROR -6


/*data structure for list reads*/
#define LIST_STRING_SIZE 80

typedef struct real_list_t{
  char name[LIST_STRING_SIZE];
  double value;
} real_list_t;

typedef struct int_list_t{
  char name[LIST_STRING_SIZE];
  int value;
} int_list_t;


int readFlashParticles(int file_identifier, double* localParticles, 
                       int localnp, int offset, int npart_props);
int readFlashLocalNP(int file_identifier, int *localnp, int *maxnp, 
		     int* start, int* end, int *tagIndex, int *npart_props,
                     int nprocs, int mype, MPI_Comm comm);
int flashFileOpen(char* filename, int* file_identifier, int create, MPI_Comm comm);
int flashFileOpenReadOnly(char* filename, hid_t* file_identifier);
void flashFileClose(int* file_identifier);


int getTimeStep(hid_t *fidIn, double *timestep);

int writeFlashParticles(int fid, double *localParticles, int localnp,
			int numPartProps, int offset);
int writeFlashParticlesHDF5(int fid, double *localParticles, int localnp,
			int numPartProps, int offset);

//int writeTrajectoryParticles(int fid, double *assigned_particles, int localnp,
//		     int numPartProps, int offset);

int trajectoryInit(hid_t fid, int numPartProps, int timesteps, int globalNumParticles); 
int copyUnknownNames(int* fidIn, int* fidOut, int myPE, int* tagIndex, int* blkIndex);
int getNumberOfParticles(hid_t *fidIn, int* numParticles, int* numPartProps);
int writeTimestepData(hid_t fidOut, int numTimeSteps, int offset, double* timesteps);
int writeTrajectoryData(hid_t *fidOut, double *particleData, int offset, int numPartProps, int timestepStart, int timeSteps, int defualtStride, int globalNumParticles, int localNumParticles);

#endif
