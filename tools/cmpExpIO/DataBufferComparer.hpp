#ifndef DATABUFFERCOMPARER_H
#define DATABUFFERCOMPARER_H

#include "mpi.h"
#include <cassert>
#include <cstdlib>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cfloat>

class DataBufferComparer {
private:
  MPI_Comm m_comm;
  int m_myPE, m_numProcs;
  mutable bool m_comparisonFailed;
  double norm(const double &a, const double &b) const;
  double max(const double &a, const double &b) const;

public:
  DataBufferComparer(const MPI_Comm &comm);
  ~DataBufferComparer();

  void DoComparison(const std::string &varName,
		    const double * const buf1,
		    const double * const buf2,
		    const size_t numDataElements) const;
  bool didComparisonFail();
};


#endif
