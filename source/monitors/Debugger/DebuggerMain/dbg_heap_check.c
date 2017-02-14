#include "dbg_heap_check.h"

/*
  This function is a wrapper around mcheck_check_all (from mcheck
  library) or dmalloc_verify (from dmalloc library) which check for
  heap corruption.  If heap corruption is detected the both functions
  abort the application.

  FLASH must be setup with the +mcheck, +dmalloc or +ddt shortcut.
*/

static int doHeapCheck = 1;

void FTOC(dbg_heap_check)(void)
{
  if (doHeapCheck) dbg_heap_check_impl();
}

void FTOC(dbg_disable_heap_check)(void)
{
  doHeapCheck = 0;
}

void FTOC(dbg_enable_heap_check)(void)
{
  doHeapCheck = 1;
}
