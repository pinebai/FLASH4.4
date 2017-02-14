#ifndef __SORTPARTICLES_H
#define __SORTPARTICLES_H

#define NONEXISTENT -1
#define DOUBLE
#ifdef DOUBLE
typedef double MYREAL;
#define MPTYPE MPI_DOUBLE_PRECISION
#else
typedef float MYREAL;
#define MPTYPE MPI_REAL
#endif
#define HUGE 99999999
#define PADDING 20

void TagListAndBounds();
int findParticleIndex(int tag);
void sortParticles(int numUnsorted);

#endif
