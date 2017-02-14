#include "WorkCoordinator.hpp"

WorkCoordinator::WorkCoordinator(int myPE, int numProcs)
{
  m_myPE = myPE;
  m_numProcs = numProcs;
}


void WorkCoordinator::GetWorkPortion(int numWorkUnits, int *workUnitStart, 
				     int *workUnitCount) const
{
  int unitBasedWorkStart, unitBasedWorkEnd;

  if ((m_numProcs > numWorkUnits) && (m_myPE == 0)) {
    std::cout << "[Info]: More processors than work units." << std::endl;
  }
  
  MPE_Decomp1d(numWorkUnits, m_numProcs, m_myPE, 
	       &unitBasedWorkStart, &unitBasedWorkEnd);

  if (unitBasedWorkStart > unitBasedWorkEnd) {
    *workUnitStart = 0;
    *workUnitCount = 0;
  } else {
    *workUnitStart = unitBasedWorkStart - 1;
    *workUnitCount = unitBasedWorkEnd - unitBasedWorkStart + 1;
  }
}


void WorkCoordinator::MPE_Decomp1d(int n, int size, int rank, 
				   int *s, int *e) const
{
  int nlocal, deficit;
  nlocal      = n / size;
  *s  = rank * nlocal + 1;
  deficit     = n % size;
  *s  = *s + ((rank < deficit) ? rank : deficit);
  if (rank < deficit) nlocal++;
  *e      = *s + nlocal - 1;
  if (*e > n || rank == size-1) *e = n;
}
