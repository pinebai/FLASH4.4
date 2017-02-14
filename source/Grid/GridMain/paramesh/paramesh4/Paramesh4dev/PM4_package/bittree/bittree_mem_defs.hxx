#ifndef _d30dc669_ee55_4dc4_b4c4_f2e8990456a9
#define _d30dc669_ee55_4dc4_b4c4_f2e8990456a9

#include "bittree_prelude.hxx"

namespace BitTree {
  struct Mem {
    size_t align, size; // size not necessarily padded to meet align

    static size_t align_size(size_t align, size_t size);
    static size_t max_align(size_t align1, size_t align2);
    static size_t min_align(size_t align1, size_t align2);
    
    static Mem zero();
    static Mem make(size_t align, size_t size);
    template<class T>
    static Mem of();
    static Mem sum(const Mem &a, const Mem &b);
    static Mem product(const Mem &a, const Mem &b);
    static Mem array(const Mem &x, size_t n);
    template<class T>
    static Mem array(size_t n);
    static Mem pad(const Mem &a);
    
    bool operator==(const Mem &that) const;
    bool operator!=(const Mem &that) const;
    
    size_t concat(const Mem &x);
  };
}
#endif
