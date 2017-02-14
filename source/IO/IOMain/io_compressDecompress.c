
#include "io_compressDecompress.h"
#include <stdlib.h>


char *io_compress(double* data, /* pointer to array of doubles*/
                    int* size,    /* how many doubles are there */
                    compress_t *comp_info /*Compression info */
                   )
{
char *comp_data;
int ctr;
double vmin,vmax,vstep;

  /* Phase1: Compute min,max,step size */
     /* Slightly complicated but more efficient min max routine */
  if (*size %2 == 0) { /* even number of elements */
     if (data[0] < data[1]) {
       vmin = data[0];
       vmax = data[1];
     } else {
       vmin = data[1];
       vmax = data[0];
     }
     ctr = 2;
  } else { /* odd number of elements */
     vmin = vmax = data[0];
     ctr = 1;
  }
  /* ctr has been set to 1 or 2 */
  for (; ctr < *size; ctr += 2) {
     if (data[ctr] < data[ctr+1]) {
        if (data[ctr] < vmin) vmin = data[ctr];
        if (data[ctr+1] > vmax) vmax = data[ctr+1];
     } else {
        if (data[ctr+1] < vmin) vmin = data[ctr+1];
        if (data[ctr] > vmax) vmax = data[ctr];
     }
  }

  comp_info->max = vmax;
  comp_info->min = vmin;
  comp_info->step = vstep = (vmax-vmin)/NUM_COMPRESSED_BLOCKS;
  comp_info->comp_size = (*size);

  /* Phase2: Allocate and populate integer array */
  comp_data = NULL;
  comp_data = (char *) malloc ( (comp_info->comp_size)* sizeof(char));

#ifdef DEBUG_IO
  if (comp_data == NULL)
     Driver_AbortC("Unable to allocate space from compressed data");
#endif

  /* Populate the array */
  for (ctr=0; ctr < *size; ctr++)
     comp_data[ctr] = ( (data[ctr]-vmin)/(vstep) + 0.5); /* Round to nearest integer */

  return comp_data;
}


/*-------------------------------------------------------------------------------*/

double *io_decompress(int* comp_data, /* pointer to array of ints*/
                    int *size, /* How many doubles in returned array */
                    compress_t *comp_info /*Compression info */
                   )
{
double *data;
double vmin,vstep;
int ctr,j,idata,base;
char cdata;

   *size = comp_info->comp_size;
   vstep = comp_info->step;
   vmin = comp_info->min;
   data = (double *) malloc( (*size) * sizeof(double));
   for (ctr=0; ctr < *size; ctr++) 
       data[ctr] = comp_data[ctr]*vstep + vmin;

   return data;
}


