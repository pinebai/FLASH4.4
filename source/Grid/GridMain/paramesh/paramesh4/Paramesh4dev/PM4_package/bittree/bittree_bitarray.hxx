#ifndef _19ce2e0d_3e4f_414d_b426_cb63ef9e4d48
#define _19ce2e0d_3e4f_414d_b426_cb63ef9e4d48

#include "bittree_bitarray_defs.hxx"
#include "bittree_bits.hxx"
#include "bittree_mem.hxx"
#include "bittree_ref.hxx"

namespace BitTree {
  template<class W>
  Ref<BitArray<W> > BitArray<W>::make(unsigned n) {
    Ref<BitArray<W> > ref;
    unsigned nw = (n + bitw)>>logw; // one extra
    Mem mem = {
      alignof(BitArray<W>),
      sizeof(BitArray<W>) + (nw-1)*sizeof(W)
    };
    BitArray<W> *it = new(ref.alloc(mem)) BitArray<W>(n);
    it->wbuf[nw-1] = W(0);
    return ref;
  }

  template<class W>
  bool BitArray<W>::get(unsigned ix) const {
    return ix < this->len ? 1 & (wbuf[ix>>logw] >> (ix & (bitw-1))) : 0;
  }

  template<class W>
  bool BitArray<W>::set(unsigned ix, bool x) {
    W w0 = wbuf[ix>>logw];
    W z = x ? W(0) : ones;
    W w1 = z^((z^w0) | one<<(ix&(bitw-1)));
    wbuf[ix>>logw] = w1;
    return w1 != w0;
  }

  template<class W>
  unsigned BitArray<W>::count() const {
    return this->count(0, this->len);
  }

  template<class W>
  unsigned BitArray<W>::count(unsigned ix0, unsigned ix1) const {
    if(ix1 <= ix0) return 0;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = (min(ix1, this->len)-1) >> logw;
    W m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((ix1-1)&(bitw-1)));
      pop += bitpop(wbuf[iw] & m);
      m = ones;
    }
    return pop;
  }

  template<class W>
  unsigned BitArray<W>::count_xor(
      const BitArray<W> *a,
      const BitArray<W> *b,
      unsigned ix0, unsigned ix1
    ) {
    if(ix1 <= ix0) return 0;
    unsigned a_ix1 = min(ix1, a->len);
    unsigned b_ix1 = min(ix1, b->len);
    unsigned a_iw1 = (a_ix1 + bitw-1) >> logw;
    unsigned b_iw1 = (b_ix1 + bitw-1) >> logw;
    unsigned iw0 = ix0 >> logw;
    unsigned iw1 = max(a_iw1, b_iw1);
    W m = ones << (ix0 & (bitw-1));
    unsigned pop = 0;
    for(unsigned iw=iw0; iw < iw1; iw++) {
      W aw = iw < a_iw1 ? a->wbuf[iw] : W(0);
      aw &= ones >> (iw+1 < a_iw1 ? 0 : bitw-1-((a_ix1-1)&(bitw-1)));
      W bw = iw < b_iw1 ? b->wbuf[iw] : W(0);
      bw &= ones >> (iw+1 < b_iw1 ? 0 : bitw-1-((b_ix1-1)&(bitw-1)));
      pop += bitpop(m & (aw ^ bw));
      m = ones;
    }
    return pop;
  }

  template<class W>
  unsigned BitArray<W>::find(unsigned ix0, unsigned nth) const {
    unsigned iw = ix0 >> logw;
    W m = ones << (ix0&(bitw-1));
    while(true) {
      if(iw >= len>>logw)
        m &= ones >> (bitw-1-((len-1)&(bitw-1)));
      W w = m & wbuf[iw];
      unsigned pop = bitpop(w);
      if(pop > nth) {
        while(nth--)
          w &= w - 1;
        return (iw<<logw) + bitffs(w)-1u;
      }
      nth -= pop;
      m = ones;
      iw += 1;
    }
  }

  template<class W>
  void BitArray<W>::fill(bool x) {
    this->fill(x, 0, this->len);
  }
  template<class W>
  void BitArray<W>::fill(bool x, unsigned ix0, unsigned ix1) {
    if(ix0 >= ix1) return;
    DBG_ASSERT(ix1 <= this->len);
    unsigned iw0 = ix0 >> logw, iw1 = (ix1-1) >> logw;
    W z = x ? W(0) : ones;
    W m = ones << (ix0 & (bitw-1));
    for(unsigned iw=iw0; iw <= iw1; iw++) {
      if(iw == iw1)
        m &= ones >> (bitw-1-((ix1-1)&(bitw-1)));
      wbuf[iw] = z ^ ((z ^ wbuf[iw]) | m);
      m = ones;
    }
  }

  template<class W>
  BitArray<W>::Reader::Reader(const BitArray<W> *a, unsigned ix0):
    a(const_cast<BitArray<W>*>(a)),
    w(a->wbuf[ix0>>logw] & (a->len < ix0+bitw ? (one<<(a->len&(bitw-1)))-1u : ~W(0))),
    ix(ix0) {
  }

  template<class W>
  template<unsigned n>
  W BitArray<W>::Reader::read() {
    unsigned n_ = n;
    BitArray<W> *a = this->a;
    W &w = this->w;
    unsigned &ix = this->ix;
    W ans;
    if(((ix+n)&(bitw-1u)) > (ix&(bitw-1u)))
      ans = (w>>(ix&(bitw-1u))) & ((one<<n_)-1u);
    else {
      ans = w>>(ix&(bitw-1u));
      w = ((ix+n)&~(bitw-1u)) < a->len ? a->wbuf[(ix>>logw)+1] : W(0);
      ans |= (w & ((one<<((ix+n)&(bitw-1u)))-1u)) << (bitw-(ix&(bitw-1u)));
    }
    ix += n;
    return ans;
  }
  
  template<class W>
  void BitArray<W>::Reader::seek(unsigned ix) {
    this->w = a->wbuf[ix>>logw] & (a->len < ix+bitw ? (one<<(a->len&(bitw-1)))-1u : ~W(0));
    this->ix = ix;
  }

  template<class W>
  BitArray<W>::Writer::Writer(BitArray<W> *host, unsigned ix0):
    BitArray<W>::Reader(host, ix0) {
  }
  template<class W>
  BitArray<W>::Writer::~Writer() {
    this->flush();
  }

  template<class W>
  template<unsigned n>
  void BitArray<W>::Writer::write(W x) {
    unsigned n_ = n;
    BitArray<W> *a = this->a;
    W &w = this->w;
    unsigned &ix = this->ix;
    DBG_ASSERT(ix + n <= a->length());
    if(((ix+n)&(bitw-1u)) > (ix&(bitw-1u))) {
      W m = ((one<<n_)-one)<<(ix&(bitw-1u));
      w = (w & ~m) | x<<(ix&(bitw-1u));
    }
    else {
      W m = ones<<(ix&(bitw-1u));
      w = (w & ~m) | x<<(ix&(bitw-1u));
      a->wbuf[ix>>logw] = w;
      w = a->wbuf[(ix>>logw)+1];
      m = (one<<((ix+n)&(bitw-1u)))-one;
      w = (w & ~m) | x>>(bitw-(ix&(bitw-1u)));
    }
    ix += n;
  }

  template<class W>
  void BitArray<W>::Writer::flush() {
    this->a->wbuf[this->ix>>logw] = this->w;
  }

  template<class W>
  FastBitArray<W>::FastBitArray(unsigned len):
    bitsref(BitArray<W>::make(len)),
    bits(bitsref),
    chksref(Ref_::new_array<unsigned>(len>>logc)),
    chks(chksref) {
  }

  template<class W>
  FastBitArray<W>::Builder::Builder(unsigned len):
    ref(Ref_::new1<FastBitArray<W>,unsigned>(len)),
    w(typename BitArray<W>::Writer(ref->bits, 0)),
    pchk(ref->chks),
    chkpop(0) {
  }

  template<class W>
  template<unsigned n>
  void FastBitArray<W>::Builder::write(W x) {
    unsigned ix = w.index();
    if((ix&(bitc-1u))+n >= bitc) {
      if(n == 1)
        *pchk++ = chkpop + x;
      else {
        W m = ~W(0) >> (BitArray<W>::bitw-(bitc-(ix&(bitc-1u))));
        *pchk++ = chkpop + bitpop(x & m);
      }
    }
    chkpop += n == 1 ? x : bitpop(x);
    w.template write<n>(x);
  }

  template<class W>
  Ref<FastBitArray<W> > FastBitArray<W>::Builder::finish() {
    w.flush();
    return ref;
  }

  template<class W>
  unsigned FastBitArray<W>::count(unsigned ix0, unsigned ix1) const {
    if(ix1>>logc > ix0>>logc) {
      unsigned pop0 = ix0 == 0 ? 0 : chks[((ix0+bitc-1u)>>logc)-1];
      unsigned pop1 = chks[(ix1>>logc)-1];
      return
        bits->count(ix0, (ix0+bitc-1u) & ~(bitc-1u)) +
        (pop1 - pop0) +
        bits->count(ix1 & ~(bitc-1u), ix1);
    }
    else
      return bits->count(ix0, ix1);
  }

  template<class W>
  unsigned FastBitArray<W>::find(unsigned ix0, unsigned nth) const {
    unsigned pop0 = this->count(0, ix0);
    unsigned a = (ix0+bitc-1u)>>logc;
    unsigned b = bits->length()>>logc;
    if(a <= b && (a==0 ? 0 : chks[a-1]) - pop0 <= nth) {
      // binary search of interval [a,b] (both inclusive)
      while(a + 10 < b) {
        unsigned x = (a+b)>>1;
        if(pop0 + nth < chks[x-1]) // x cant be 0
          b = x-1;
        else
          a = x;
      }
      // linear search
      while(a < b && chks[a]-pop0 <= nth)
        a += 1;
      nth -= (a==0 ? 0 : chks[a-1]) - pop0;
      return bits->find(a<<logc, nth);
    }
    else
      return bits->find(ix0, nth);
  }
}
#endif
