#include "chombo_interface.hpp"

void Chombo_Interface::InitGrid(size_t v, size_t x, size_t y, size_t z)
{
  size_t i;
  m_pMem = (double*) malloc(v * x * y * z * sizeof(double));
  for (i=0; i<(v*x*y*z); ++i) {
    m_pMem[i] = (double) i;
  }
}

double * Chombo_Interface::GetDataPtr(size_t blkID)
{
  if (blkID == 1) {
    return m_pMem;
  } else {
    return NULL;
  }
}

void Chombo_Interface::FreeGrid()
{
  free (m_pMem);
}
