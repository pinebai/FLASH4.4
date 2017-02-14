#ifndef CHOMBO_INTERFACE_HPP
#define CHOMBO_INTERFACE_HPP
#include <cstdlib>

class Chombo_Interface {
private:
  double *m_pMem;
public:
  void InitGrid(size_t v, size_t x, size_t y, size_t z);
  double * GetDataPtr(size_t blkID);
  void FreeGrid();
};

#endif
