#ifndef WORKCOORDINATOR_H
#define WORKCOORDINATOR_H

#include <iostream>
#include <cstdlib>
#include "mpi.h"

class WorkCoordinator {
private:
  void MPE_Decomp1d(int n, int size, int rank, int *s, int *e) const;
  int m_myPE;
  int m_numProcs;

public:
  WorkCoordinator(int myPE, int numProcs);
  void GetWorkPortion(int numWorkUnits, int *workUnitStart, int *workUnitCount) const;
};

#endif
