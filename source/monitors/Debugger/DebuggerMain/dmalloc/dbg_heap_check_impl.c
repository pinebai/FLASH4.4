#include "dbg_heap_check_impl.h"

/*
  This function is a wrapper around dmalloc_verify which checks for
  heap corruption.  If heap corruption is detected the dmalloc library
  aborts.

  FLASH must be setup with either the +dmalloc or +ddt shortcut.
*/

void dbg_heap_check_impl(void)
{
#ifdef DMALLOC
  int i;
  i = dmalloc_verify(NULL);
#endif
}
