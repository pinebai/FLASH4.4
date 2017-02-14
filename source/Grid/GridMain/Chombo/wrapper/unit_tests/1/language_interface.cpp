#include <cstdlib>
#include "chombo_interface.hpp"
#include "Flash.h"

/* The functions in this file operate on our Chombo interface
   object.  In this object we actually interact with Chombo library */
Chombo_Interface *pChomboInterface;


extern "C" void language_interface_init()
{
  pChomboInterface = new Chombo_Interface();
  pChomboInterface->InitGrid(NVAR,NXB,NYB,NZB);
}


extern "C" void language_interface_get_blk_ptr(int blkID, void **ppv)
{
  *ppv = (void*) pChomboInterface->GetDataPtr((size_t)blkID);
}


extern "C" void language_interface_finalize()
{
  pChomboInterface->FreeGrid();
  delete(pChomboInterface);
}
