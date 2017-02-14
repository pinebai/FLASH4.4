#ifndef _f098c2f1_fd19_4bfe_93bc_2f8020441d84
#define _f098c2f1_fd19_4bfe_93bc_2f8020441d84

#include "bittree_mem_defs.hxx"

#include <algorithm>

namespace BitTree {
  using namespace std;
  
  inline size_t Mem::align_size(size_t align, size_t size) {
    return (size + (align-1)) & ~(align-1);
  }
  inline size_t Mem::max_align(size_t align1, size_t align2) {
    return ((align1-1) | (align2-1)) + 1;
  }
  inline size_t Mem::min_align(size_t align1, size_t align2) {
    return ((align1-1) & (align2-1)) + 1;
  }
  inline Mem Mem::make(size_t align, size_t size) {
    Mem m = { align, size };
    return m;
  }
  inline Mem Mem::zero() {
    Mem m = { 1, 0 };
    return m;
  }
  template<class T>
  inline Mem Mem::of() {
    Mem m = { alignof(T), sizeof(T) };
    return m;
  }
  inline Mem Mem::pad(const Mem &x) {
    Mem m = { x.align, align_size(x.align, x.size) };
    return m;
  }
  inline Mem Mem::sum(const Mem &a, const Mem &b) {
    size_t align = max_align(a.align, b.align);
    Mem m = { align, align_size(align, max(a.size, b.size)) };
    return m;
  }
  inline Mem Mem::product(const Mem &a, const Mem &b) {
    Mem m = { max_align(a.align, b.align), align_size(b.align, a.size) + b.size };
    return m;
  }
  inline Mem Mem::array(const Mem &x, size_t n) {
    Mem m = { n ? x.align : 1, n*align_size(x.align, x.size) };
    return m;
  }
  template<class T>
  inline Mem Mem::array(size_t n) {
    return Mem::array(Mem::of<T>(), n);
  }
  inline bool Mem::operator==(const Mem &that) const {
    return this->align == that.align && this->size == that.size;
  }
  inline bool Mem::operator!=(const Mem &that) const {
    return this->align != that.align || this->size != that.size;
  }
  inline size_t Mem::concat(const Mem &x) {
    size_t p = align_size(x.align, size);
    align = max_align(align, x.align);
    size = p + x.size;
    return p;
  }	
}
#endif
