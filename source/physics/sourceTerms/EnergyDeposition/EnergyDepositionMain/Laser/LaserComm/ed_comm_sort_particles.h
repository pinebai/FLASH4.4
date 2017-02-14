#ifndef ED_COMM_SORT_PARTICLES_H
#define ED_COMM_SORT_PARTICLES_H

#include <stdio.h>
#include <stdlib.h>

#ifdef UNIT_TEST
# define RAY_ATTR_COUNT 3
# define RAY_POSX 1
# define RAY_POSY 2
# define RAY_BLCK 3
# define RAY_COUNT 8
# define FTOC(x) x
#else
# include "EnergyDeposition.h"
# include "mangle_names.h"
#endif

int ascending_block_order(const void *,
			  const void *);
void FTOC(ed_comm_sort_particles)(double [],
				  const int * const);
void print_ray_blockids(const double rays[],
			const int rayCount);

#endif
