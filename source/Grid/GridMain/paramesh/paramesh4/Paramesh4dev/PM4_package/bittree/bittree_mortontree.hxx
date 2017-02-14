#ifndef _328c0dfe_e65a_477f_b38d_69b42f91736f
#define _328c0dfe_e65a_477f_b38d_69b42f91736f

#include "bittree_mortontree_defs.hxx"
#include "bittree_bitarray.hxx"
#include "bittree_bits.hxx"
#include "bittree_mem.hxx"
#include "bittree_ref.hxx"

#include <algorithm>
#include <limits>
#include <iomanip>

namespace BitTree {
  template<unsigned D>
  unsigned rect_coord_to_mort(const unsigned domain[D], const unsigned coord[D]) {
    unsigned x[D], box[D];
    for(unsigned d=0; d < D; d++) {
      x[d] = coord[d];
      box[d] = domain[d];
    }
    unsigned mort = 0u;
    // bisect box until it's just one element
    while(true) {
      // find dim that can fit biggest pow2 strictly inside box
      unsigned max_pow2 = 0u;
      unsigned max_d;
      for(unsigned d=0; d < D; d++) {
        unsigned p2 = glb_pow2(box[d]-1u);
        if(p2 >= max_pow2) {
          max_pow2 = p2;
          max_d = d;
        }
      }
      if(max_pow2 == 0u)
        return mort; // the box is just one, we're done
      if(x[max_d] < max_pow2)
        box[max_d] = max_pow2;
      else {
        unsigned pop = 1u;
        for(unsigned d=0; d < D; d++)
          pop *= d == max_d ? max_pow2 : box[d];
        mort += pop;
        x[max_d] -= max_pow2;
        box[max_d] -= max_pow2;
      }
    }
  }
  
  template<unsigned D>
  void rect_mort_to_coord(const unsigned domain[D], unsigned mort, unsigned coord[D]) {
    unsigned box[D];
    for(unsigned d=0; d < D; d++) {
      coord[d] = 0u;
      box[d] = domain[d];
    }
    // bisect box until it's just one element
    while(true) {
      // find dim that can fit biggest pow2 strictly inside box
      unsigned max_pow2 = 0u;
      unsigned max_d;
      for(unsigned d=0; d < D; d++) {
        unsigned p2 = glb_pow2(box[d]-1u);
        if(p2 >= max_pow2) {
          max_pow2 = p2;
          max_d = d;
        }
      }
      if(max_pow2 == 0u)
        return; // the box is just one, we're done
      unsigned pop = 1u; // left sub-box population
      for(unsigned d=0; d < D; d++)
        pop *= d == max_d ? max_pow2 : box[d];
      if(mort < pop)
        box[max_d] = max_pow2;
      else {
        mort -= pop;
        coord[max_d] += max_pow2;
        box[max_d] -= max_pow2;
      }
    }
  }
  
  template<unsigned D, class W>
  Ref<MortonTree<D,W> > MortonTree<D,W>::make(
      const unsigned size[D],
      const bool includes[]
    ) {
    Ref<MortonTree<D,W> > ref;
    Mem mem = {
      alignof(MortonTree<D,W>),
      sizeof(MortonTree<D,W>) + 0*sizeof(Level)
    };
    MortonTree<D,W> *it = new(ref.alloc(mem)) MortonTree<D,W>();
    
    unsigned blkpop = 1;
    for(unsigned d=0; d < D; d++) {
      it->lev0_blks[d] = size[d];
      blkpop *= size[d];
    }
    it->levs = 1;
    // the first blkpop bits of our bitarray are 'inclusion' bits
    // after that come the actual block bits for all but the last level, which has no bits
    it->id0 = blkpop;
    unsigned lev0_id1 = it->id0;
    typename FastBitArray<W>::Builder bldr(blkpop);
    // generate inclusion bits
    for(unsigned mort=0; mort < blkpop; mort++) {
      unsigned x[D];
      rect_mort_to_coord<D>(size, mort, x);
      unsigned ix = 0;
      for(unsigned d=D; d--;)
        ix = size[d]*ix + x[d];
      unsigned include = includes[ix] ? 1 : 0;
      lev0_id1 += include;
      bldr.template write<1>(include);
    }
    // ok bitarray done.  since there's only one level, we dont store any block bits
    it->bits = bldr.finish();
    it->level[0].id1 = lev0_id1;
    return ref;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::levels() const {
    return this->levs;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::blocks() const {
    return this->level[this->levs-1].id1 - this->id0;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::top_size(unsigned dim) const {
    return this->lev0_blks[dim];
  }
  
  template<unsigned D, class W>
  unsigned MortonTree<D,W>::id_upper_bound() const {
    return this->level[this->levs-1].id1;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::level_id0(unsigned lev) const {
    return lev == 0 ? this->id0 : this->level[lev-1].id1;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::level_blocks(unsigned lev) const {
    return this->level[lev].id1 - (lev == 0 ? this->id0 : this->level[lev-1].id1);
  }

  template<unsigned D, class W>
  bool MortonTree<D,W>::block_is_parent(unsigned id) const {
    return levs > 1 && id < level[levs-2].id1 && bits->get(id);
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::block_level(unsigned id) const {
    unsigned lev = 0;
    while(level[lev].id1 <= id)
      lev += 1;
    return lev;
  }

  template<unsigned D, class W>
  template<class X>
  bool MortonTree<D,W>::inside(unsigned lev, const X x[D]) const {
    unsigned x0[D];
    for(unsigned d=0; d < D; d++) {
      x0[d] = unsigned(x[d] >> lev);
      if(x0[d] >= lev0_blks[d])
        return false;
    }
    return bits->get(rect_coord_to_mort<D>(lev0_blks, x0));
  }
  
  template<unsigned D, class W>
  template<class X>
  typename MortonTree<D,W>::template Block<X>
  MortonTree<D,W>::identify(unsigned lev, const X x[D]) const {
    const unsigned levs = this->levs;
    const unsigned id0 = this->id0;
    const FastBitArray<W> *fast_bits = this->bits; // use this for counting
    const BitArray<W> *bits = fast_bits->bit_array(); // use this for bit access
    Block<X> ans;
    unsigned ix; // index of current block in current level
    { // top level=0
      unsigned x0[D];
      for(unsigned d=0; d < D; d++) {
        x0[d] = unsigned(x[d] >> lev); // coarsen x to top level
        DBG_ASSERT(x0[d] < lev0_blks[d]);
        ans.coord[d] = X(x0[d]);
      }
      ix = rect_coord_to_mort<D>(lev0_blks, x0);
      DBG_ASSERT(bits->get(ix));
      ix = fast_bits->count(0, ix); // discount excluded blocks
    }
    ans.mort = 0;
    // bisection iteration
    for(unsigned a_lev=0; a_lev < levs; a_lev++) {
      unsigned a_id0 = a_lev==0 ? id0 : level[a_lev-1].id1;
      unsigned a_id = a_id0 + ix;
      ans.mort += ix;
      unsigned inside = 0u;
      bool is_par = a_lev+1u < levs && bits->get(a_id);
      if(is_par && a_lev < lev) { // if we're bisecting further
        for(unsigned d=0; d < D; d++) {
          X xd = x[d] >> (lev-a_lev-1u);
          ans.coord[d] <<= 1;
          if(xd >= ans.coord[d]+1u) {
            ans.coord[d] += 1u;
            inside += 1u << d;
            ans.mort += d == D-1 ? 1 : 0;
          }
        }
      }
      else if(a_lev <= lev) { // we have the result block
        lev = a_lev; // stop this from running again
        ans.id = a_id;
        ans.level = a_lev;
        ans.is_parent = is_par;
        if(is_par) inside = 1u<<(D-1); // include first half of children
      }
      unsigned parbef = this->parents_before(a_lev, ix);
      ix = (parbef<<D) + inside;
    }
    return ans;
  }

  template<unsigned D, class W>
  template<class X>
  typename MortonTree<D,W>::template Block<X>
  MortonTree<D,W>::locate(unsigned id) const {
    const unsigned id0 = this->id0;
    Block<X> ans;
    ans.id = id;
    ans.mort = 0;
    ans.is_parent = this->block_is_parent(id);
    unsigned lev = 0;
    while(level[lev].id1 <= id)
      lev += 1;
    ans.level = lev;
    for(unsigned d=0; d < D; d++)
      ans.coord[d] = X(0u);
    // index on this level
    unsigned ix = id - (lev == 0 ? id0 : level[lev-1].id1);
    { // count children of all preceeding parents in morton index
      unsigned down = ix;
      for(unsigned lev1=lev; lev1 < levs; lev1++) {
        down = this->parents_before(lev1, down) << D;
        if(lev1 == lev && ans.is_parent)
          down += 1u<<(D-1);
        ans.mort += down;
      }
    }
    // walk up the levels
    while(0 < lev) {
      for(unsigned d=0; d < D; d++)
        ans.coord[d] += X(ix>>d & 1u) << (ans.level-lev);
      ans.mort += ix + (ix>>(D-1) & 1u);
      ix = this->parent_find(lev-1, ix>>D) - (lev-1==0 ? id0 : level[lev-2].id1);
      lev -= 1;
    }
    ans.mort += ix;
    { // top level=0
      unsigned x0[D];
      ix = this->bits->find(0, ix); // account for excluded blocks
      rect_mort_to_coord<D>(lev0_blks, ix, x0);
      for(unsigned d=0; d < D; d++)
        ans.coord[d] += X(x0[d]) << ans.level;
    }
    return ans;
  }

  template<unsigned D, class W>
  Ref<MortonTree<D,W> > MortonTree<D,W>::refine(
      Ref<BitArray<W> > delta_
    ) const {
    using namespace std;
    
    const unsigned a_levs = this->levs;
    const BitArray<W> *a_bits = this->bits->bit_array();
    const unsigned id0 = this->id0;
    const BitArray<W> *delta = delta_;
    
    // count the new number of levels, blocks, and bits
    unsigned b_id1 = this->level[0].id1;
    unsigned b_bitlen = id0;
    unsigned b_levs = 1;
    for(unsigned lev=0; lev < a_levs; lev++) {
      unsigned lev_id0 = lev == 0 ? id0 : this->level[lev-1].id1;
      unsigned lev_id1 = this->level[lev].id1;
      unsigned b_pars = BitArray<W>::count_xor(a_bits, delta, lev_id0, lev_id1);
      if(b_pars != 0) b_bitlen = b_id1;
      b_id1 += b_pars << D;
      if(b_pars == 0) break;
      b_levs += 1;
    }
    
    // new bit tree
    Ref<MortonTree<D,W> > b_ref_tree;
    MortonTree<D,W> *b_tree; {
      Mem m = {
        alignof(MortonTree<D,W>),
        sizeof(MortonTree<D,W>) + (b_levs-1)*sizeof(Level)
      };
      b_tree = new(b_ref_tree.alloc(m)) MortonTree<D,W>();
      b_tree->levs = b_levs;
      b_tree->id0 = id0;
      for(unsigned d=0; d < D; d++)
        b_tree->lev0_blks[d] = this->lev0_blks[d];
      // still must initialize b_tree->bits
    }
    
    // apply delta and insert/remove blocks
    typename BitArray<W>::Reader a_r(a_bits), del_r(delta, id0);
    typename FastBitArray<W>::Builder b_w(b_bitlen);
    
    // copy inclusion bits
    while(a_r.index() < id0)
      b_w.template write<1>(a_r.template read<1>());
    
    // do level 0...
    b_tree->level[0].id1 = this->level[0].id1;
    while(b_w.index() < b_bitlen && a_r.index() < this->level[0].id1) {
      // apply delta
      b_w.template write<1>(a_r.template read<1>() ^ del_r.template read<1>());
    }
    
    // readers of previous level
    typename BitArray<W>::Reader a_rp(a_bits, id0), del_rp(delta, id0);
    // do remaining levels
    unsigned lev = 1;
    while(b_w.index() < b_bitlen) {
      if(a_rp.template read<1>()) { // it was a parent
        bool still_a_parent = !del_rp.template read<1>();
        // read kids, apply delta
        W b_kids = a_r.template read<(1<<D)>() ^ del_r.template read<(1<<D)>();
        // if it became a leaf then we just dont write out the kids
        if(still_a_parent)
          b_w.template write<(1<<D)>(b_kids);
      }
      else { // it was a leaf
        if(del_rp.template read<1>()) // and it became a parent!
          b_w.template write<(1<<D)>(0);
      }
      if(a_rp.index() == this->level[lev-1].id1 || b_w.index() == b_bitlen) {
        b_tree->level[lev].id1 = b_w.index();
        lev += 1;
      }
    }
    
    b_tree->level[b_levs-1].id1 = b_id1;
    b_tree->bits = b_w.finish();
    return b_ref_tree;
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::parents_before(unsigned lev, unsigned ix) const {
    if(lev >= levs-1) return 0;
    unsigned id0 = lev == 0 ? this->id0 : level[lev-1].id1;
    return bits->count(id0, id0 + ix);
  }

  template<unsigned D, class W>
  unsigned MortonTree<D,W>::parent_find(unsigned lev, unsigned par_ix) const {
    unsigned id0 = lev == 0 ? this->id0 : level[lev-1].id1;
    return bits->find(id0, par_ix);
  }

  template<class W>
  std::ostream& operator<<(std::ostream &o, const MortonTree<2,W> *x) {
    using namespace std;
    
    unsigned levs = x->levels();
    for(unsigned lev=0; lev < levs; lev++) {
      cout << "lev=" << lev << '\n';
      unsigned ij[2];
      for(ij[1]=0; ij[1] < x->top_size(1)<<lev; ij[1]++) {
        for(ij[0]=0; ij[0] < x->top_size(0)<<lev; ij[0]++) {
          if(x->inside(lev, ij)) {
            typename MortonTree<2,W>::template Block<unsigned> b0 = x->template identify<unsigned>(lev, ij);
            typename MortonTree<2,W>::template Block<unsigned> b1 = x->template locate<unsigned>(b0.id);
            DBG_ASSERT(b0.id == b1.id);
            DBG_ASSERT(b0.level == b1.level);
            DBG_ASSERT(b0.mort == b1.mort);
            DBG_ASSERT(b0.coord[0] == b1.coord[0] && b0.coord[1] == b1.coord[1]);
            o << std::setw(4) << b0.mort;
          }
          else
            o << std::setw(4) << ' ';
        }
        o << '\n';
      }
    }
    return o;
  }
}
#endif
