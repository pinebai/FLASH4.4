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

int readFlashParticles(int file_identifier, double* localParticles, 
                       int localnp, int offset, int npart_props);
int readFlashLocalNP(int file_identifier, int *localnp, int *maxnp, 
		     int* start, int* end, int *tagIndex, int *npart_props,
                     int nprocs, int mype, MPI_Comm comm);
int flashFileOpen(char* filename, int* file_identifier);
void flashFileClose(int* file_identifier);

int writeFlashParticles(int fid, double *localParticles, int localnp,
			int numPartProps, int offset);
int writeFlashParticlesHDF5(int fid, double *localParticles, int localnp,
			int numPartProps, int offset);


#endif
