#include <math.h>
#include <stdio.h>

#include "flash_reader.h"
#include "sameblock.h"
#include "options.h"

#define SQUARE(X) (X)*(X)

int sameblock(int dim, options_t* opts, FR_File *A, FR_File *B, 
	      int NA, int NB){
  
  int d, a, b;
  double dist;

  if(opts->ignore==0){ 
    if(A->lref[NA]!=B->lref[NB])
      return 0;
  }

  if(opts->mesh_tol==0.){ /*bit-level check*/
   
    /* check bounding box */
    for(d=0; d<2*dim; d++){
      a=2*dim*NA + d;
      b=2*dim*NB + d;
      if(A->bbox[a]!=B->bbox[b])
	return 0;
    }

    /*check coord and size (redundant in theory, but...) */
    for(d=0; d<dim; d++){
      a=dim*NA+d;
      b=dim*NB+d;
      
      if((A->coord[a]!=B->coord[b]) || \
	 (A->size[a]!=B->size[b]))
	return 0;
      
    }
    return 1;
  }

  /* Actually use opts->mesh_tol */
  dist = 0;
  for(d=0; d<dim; d++){
    a=dim*NA+d;
    b=dim*NB+d;
    dist += SQUARE(A->coord[a] - B->coord[b]) / \
            SQUARE(A->size[a] + B->size[b]);
  }
  dist = sqrt(dist);
  if(opts->verbose){
    printf("Norm between blocks %d and %d: %g ", NA, NB, dist);
    if(dist < opts->mesh_tol)
      printf("(identified)\n");
    else
      printf("(different)\n");
  }
  return (dist < opts->mesh_tol);
}
