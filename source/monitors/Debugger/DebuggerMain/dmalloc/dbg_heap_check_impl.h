#include <stdlib.h>
#include "mangle_names.h"

#ifdef DMALLOC
# include "dmalloc.h"
#else

# ifdef DDT
/* Under the hood DDT is using the dmalloc library, but does not
   expose any of the internals to us.  I copy the relevant code from
   dmalloc.h (this is an ugly hack) */
#  define DMALLOC

/* this defines what type the standard void memory-pointer is */
#  if (defined(__STDC__) && __STDC__ == 1) || defined(__cplusplus) || defined(STDC_HEADERS) || defined(_ISO_STDLIB_ISO_H)
#   define DMALLOC_PNT             void *
#  else
#   define DMALLOC_PNT             char *
#  endif

extern int dmalloc_verify(const DMALLOC_PNT pnt);

# endif
#endif

void dbg_heap_check_impl(void);
