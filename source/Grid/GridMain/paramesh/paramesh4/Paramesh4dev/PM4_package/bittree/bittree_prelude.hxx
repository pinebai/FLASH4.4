#ifndef _36e30480_1361_44b3_8616_ab7c8bae0c41
#define _36e30480_1361_44b3_8616_ab7c8bae0c41

#include <cstddef>
#include <cstdlib>
#include <climits>
#include <exception>
#include <cstring>
#include <iostream>

#define alignof __alignof__
#define restrict __restrict__

#define DBG_ASSERT(x) if(!(x)) { std::cout << "FAILED " __FILE__ "@" << __LINE__ << ": " #x << std::endl; abort(); }
#define RT_ASSERT(x) if(!(x)) { std::cout << "FAILED " __FILE__ "@" << __LINE__ << ": " #x << std::endl; abort(); }

namespace BitTree {
  template<size_t b, size_t x> struct Log { enum { val = 1 + Log<b,x/b>::val }; };
  template<size_t b> struct Log<b,1> { enum { val = 0 }; };
  template<size_t b> struct Log<b,0> { enum { val = 0 }; };
}
#endif
