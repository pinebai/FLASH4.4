
/*********************

io_compress takes an array of doubles and returns an array of chars (via return value)
   and other information in compress_info (modifiable argument)

io_decompress takes an array of chars and compress_info and 
   returns an array of doubles (via return value)

********************/



/* 
   This function takes an array of doubles 
   and returns an array of chars of same size 
   together the compression info
*/

#define NUM_COMPRESSED_BLOCKS (256)

typedef struct compress_t {
  double min,max,step;
  int comp_size,local_blocks;
} compress_t;

char *io_compress(double* data, /* pointer to array of doubles*/
                    int* size,    /* how many doubles are there */
                    compress_t *comp_info /*Compression info */
                   );

/* 
   This function takes an array of integers and scaling parameters
   and returns an array of doubles (of 4 times the size)
*/

double *io_decompress(int* comp_data, /* pointer to array of ints*/
                    int *size, /* How many doubles in returned array */
                    compress_t *comp_info /*Compression info */
                   );

