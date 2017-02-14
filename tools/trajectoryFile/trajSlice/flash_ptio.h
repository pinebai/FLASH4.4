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


int fileOpenRead(char *filename, hid_t *fid);
int fileOpenWrite(char *filename, hid_t *fid);

int fileClose(hid_t *fid);

int copyUnknownNames(hid_t fidIn, hid_t fidOut, int* tagIndex, int* blkIndex,
		     int useThreshold, char *thresholdName, 
		     int *thresholdIndex);

int copyTimesteps(hid_t fidIn, hid_t fidOut);
int copyParticleSelection(hid_t fileIn, hid_t fileOut, 
			  int partStart, int partStride, int partCount);
int copyParticleThreshold(hid_t fileIn, hid_t fileOut, 
			  int maxCount, int maxThreshold, int thresholdIndex, 
			  double thresholdValue);
int findThresholdTrajectories(hid_t fid, int thresholdIndex, 
			      double thresholdValue, 
			      int *partCount, int(*compare)(double, double));
int compareLess(double a, double b);
int compareGreater(double a, double b);



#endif
