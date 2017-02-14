/* String comparison functions, like str[n]cmp, but for names consisting of one word which may be terminated by
   either a SPACE or a NULL character.
These return 0 for strings considered matching, -1 for strings considered different.
*/

#include <string.h>
#include "namecmp.h"

int namecmp(const char *s1, const char *s2){
  
  for ( ; (*s1) && (*s1 != ' '); s1++, s2++) {
    if (*s1 != *s2) {
      return -1;
    }
  }
  /* if s2 has not terminated, it's not a match */
  if (*s2 && *s2 != ' ') return -1;
  else return 0;
}

int namencmp(const char *s1, const char *s2, size_t n){
  
  size_t i = 0;
  for ( ; (i < n) && (*s1) && (*s1 != ' '); s1++, s2++, i++) {
    if (*s1 != *s2) {
      return -1;
    }
  }
  if (i==n) return 0;
  /* if s2 has not terminated, it's not a match */
  if (*s2 && *s2 != ' ') return -1;
  else return 0;
}
