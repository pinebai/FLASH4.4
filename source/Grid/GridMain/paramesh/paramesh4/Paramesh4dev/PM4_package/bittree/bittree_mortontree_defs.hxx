#ifndef _9aa50c52_f364_47ed_876f_ef6dee01de72
#define _9aa50c52_f364_47ed_876f_ef6dee01de72

#include "bittree_prelude.hxx"
#include "bittree_bitarray_defs.hxx"
#include "bittree_ref_defs.hxx"

#include <iostream>

namespace BitTree {
  template<unsigned D>
  unsigned rect_coord_to_mort(const unsigned domain[D], const unsigned coord[D]);
  
  template<unsigned D>
  void rect_mort_to_coord(const unsigned domain[D], unsigned mort, unsigned coord[D]);
  
  template<unsigned D, class W>
  class MortonTree {
    unsigned levs;
    unsigned lev0_blks[D];
    Ref<FastBitArray<W> > bits;
    unsigned id0; // id of first block
    struct Level {
      unsigned id1; // exclusive upper bound on block ids for this level
    } level[1];/*[levs]*/
  private:
    MortonTree() {}
  public:
    template<class X>
    struct Block {
      unsigned id;
      unsigned mort;
      unsigned level;
      bool is_parent;
      X coord[D];
    };
  public:
    static Ref<MortonTree<D,W> > make(
      const unsigned blks[D],
      const bool includes[]=0x0
    );
    unsigned levels() const;
    unsigned blocks() const;
    unsigned top_size(unsigned dim) const;
    unsigned id_upper_bound() const;
    unsigned level_id0(unsigned lev) const;
    unsigned level_blocks(unsigned lev) const;
    
    bool block_is_parent(unsigned id) const;
    unsigned block_level(unsigned id) const;
    
    template<class X>
    bool inside(unsigned lev, const X coord[D]) const;
    template<class X>
    Block<X> identify(unsigned lev, const X coord[D]) const;
    template<class X>
    Block<X> locate(unsigned id) const;
    
    Ref<MortonTree<D,W> > refine(Ref<BitArray<W> > delta) const;
  private:
    unsigned parents_before(unsigned lev, unsigned ix) const;
    unsigned parent_find(unsigned lev, unsigned par_ix) const;
  };

  template<class W>
  std::ostream& operator<<(std::ostream &o, const MortonTree<2,W> *x);
}
#endif
