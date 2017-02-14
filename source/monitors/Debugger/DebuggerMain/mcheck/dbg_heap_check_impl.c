#include "dbg_heap_check_impl.h"

/*
  This function is a wrapper around mcheck_check_all which checks for
  heap corruption.  If heap corruption is detected the mcheck library
  aborts.

  FLASH must be setup with the +mcheck shortcut.
*/

void dbg_heap_check_impl(void)
{
#ifdef MCHECK
  mcheck_check_all();
#endif
}
