#ifndef _a3ff44bb_f654_43ab_9ed5_642f0748b709
#define _a3ff44bb_f654_43ab_9ed5_642f0748b709

#include "bittree_prelude.hxx"

#if defined(__IBMCPP__)
# include "builtins.h"
#endif

namespace BitTree {
#if defined(__IBMCPP__)
  inline int bitpop(unsigned int x) {
    return __popcnt4(x);
  }
  inline int bitpop(unsigned long x) {
    return __popcnt4(x);
  }
  inline int bitpop(unsigned long long x) {
    return __popcnt8(x);
  }
#elif defined(__GNUC__)
  inline int bitpop(unsigned int x) {
    return __builtin_popcount(x);
  }
  inline int bitpop(unsigned long x) {
    return __builtin_popcountl(x);
  }
  inline int bitpop(unsigned long long x) {
    return __builtin_popcountll(x);
  }
#else
  template<class X> // X should be unsigned
  inline int bitpop(X x) {
    for(unsigned i=0; i < Log<2,sizeof(X)*CHAR_BIT>::val; i++) {
      unsigned b = 1u<<i;
      X m = (~X(0u))/((1u<<b)+1u);
      x = (x & m) + (x>>b & m);
    }
    return int(x);
  }
#endif

#if defined(__IBMCPP__)
  inline int bitffs(unsigned int x) {
    return x == 0 ? 0 : 1 + __cnttz4(x);
  }
  inline int bitffs(unsigned long x) {
    return x == 0 ? 0 : 1 + __cnttz4(x);
  }
  inline int bitffs(unsigned long long x) {
    return x == 0 ? 0 : 1 + __cnttz8(x);
  }
#elif defined(__GNUC__)
  inline int bitffs(unsigned int x) {
    return __builtin_ffs(x);
  }
  inline int bitffs(unsigned long x) {
    return __builtin_ffsl(x);
  }
  inline int bitffs(unsigned long long x) {
    return __builtin_ffsll(x);
  }
#else
  template<class X> // X should be unsigned
  inline int bitffs(X x) {
    int b = 0;
    for(int i=Log<2,sizeof(X)*CHAR_BIT>::val; i--;) {
      int s = x & ((X(1u)<<(1u<<i))-1u) ? 0 : 1<<i;
      b += s;
      x >>= s;
    }
    return b + (x ? 1 : 0);
  }
#endif

  // returns the greatest power of 2 less-or-equal to x, and 0 if x=0
  template<class X> // X should be unsigned
  inline X glb_pow2(X x) {
    // should unroll
    for(unsigned i=0; i < Log<2,sizeof(X)*CHAR_BIT>::val; i++)
      x |= x >> (1<<i);
    return x - (x>>1);
  }

  // returns the least power of 2 greater-or-equal to x,
  // and 0 if x=0 or if the answer exceeds the numerical range of type X
  template<class X> // X should be unsigned
  inline X lub_pow2(X x) {
    x -= 1;
    // should unroll
    for(unsigned i=0; i < Log<2,sizeof(X)*CHAR_BIT>::val; i++)
      x |= x >> (1<<i);
    return x + 1;
  }
}
#endif
