#include "ed_comm_sort_particles.h"

void FTOC(ed_comm_sort_particles)(double rays[], const int * const pRayCount)
{
  const int rayCount = *pRayCount;
  if (rayCount >= 0) {
    qsort(rays,
          rayCount,
          RAY_ATTR_COUNT*sizeof(double),
          ascending_block_order);
  }
}

int ascending_block_order(const void * v1, const void * v2)
{
  const double blk1 = ((const double*)v1)[RAY_BLCK-1];
  const double blk2 = ((const double*)v2)[RAY_BLCK-1];

  if (blk1 >= 1.0) {
    if (blk2 >= 1.0) {
      /* blk1 and blk2 are valid block IDs */
      if (blk1 < blk2) return -1;
      if (blk1 > blk2) return 1;
      if (blk1 == blk2) return 0;
    } else {
      /* blk1 is the only valid block ID - blk1 goes first */
      return -1;
    }
  } else {
    if (blk2 >= 1.0) {
      /* blk2 is the only valid block ID - blk2 goes first */
      return 1;
    } else {
      /* blk1 and blk2 are invalid block IDs */
      if (blk1 < blk2) return 1;
      if (blk1 > blk2) return -1;
      if (blk1 == blk2) return 0;
    }
  }

  fprintf(stderr, "Sorting fail!!! blk1=%.0lf, blk2=%.0lf.\n", blk1, blk2);
  exit(EXIT_FAILURE);
}

void print_ray_blockids(const double rays[], const int rayCount)
{
   int i;
   for (i=0; i<rayCount; i++)
     printf("Ray %d:  Block %.0lf\n",
            i+1, rays[(i*RAY_ATTR_COUNT)+(RAY_BLCK-1)]);
}


#ifdef UNIT_TEST
int main()
{
  double ed_rays[RAY_ATTR_COUNT*RAY_COUNT] = {
    1.0, 1.0, 3.0,
    2.0, 2.0, 1.0,
    3.0, 3.0, -1.0,
    4.0, 4.0, 4.0,
    5.0, 5.0, 2.0,
    6.0, 6.0, -1.0,
    7.0, 7.0, -12345.0,
    8.0, 8.0, 3.0
  };
  const int rayCount = RAY_COUNT;

  printf("List before sorting:\n");
  print_ray_blockids(ed_rays, rayCount);

  FTOC(ed_comm_sort_particles)(ed_rays, &rayCount);

  printf("\nList after sorting:\n");
  print_ray_blockids(ed_rays, rayCount);

  return 0;
}
#endif
