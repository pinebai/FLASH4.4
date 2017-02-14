#ifndef _1345010d_9317_43fa_90c7_d2dc7bc16926
#define _1345010d_9317_43fa_90c7_d2dc7bc16926

#include "bittree_prelude.hxx"
#include "bittree_ref_defs.hxx"

namespace BitTree {
  template<class W>
  class BitArray {
  public:
    static const unsigned logw = Log<2,sizeof(W)*CHAR_BIT>::val;
    static const unsigned bitw = 1u << logw;
    static const W one = W(1), ones = ~W(0);
  private:
    unsigned len;
    W wbuf[1];
  private:
    BitArray(unsigned len): len(len) {}
  public:
    static Ref<BitArray<W> > make(unsigned n);
  public:
    unsigned length() const { return len; }
    unsigned word_count() const { return (len+bitw-1u)>>logw; }
    W* word_buf() { return wbuf; }
    
    bool get(unsigned ix) const;
    bool set(unsigned ix, bool x);
    
    // count 1's in interval [ix0,ix1)
    unsigned count(unsigned ix0, unsigned ix1) const;
    // count 1's in whole array
    unsigned count() const;
    
    static unsigned count_xor(
      const BitArray<W> *a,
      const BitArray<W> *b,
      unsigned ix0, unsigned ix1
    );
    
    unsigned find(unsigned ix0, unsigned nth) const;
    
    void fill(bool x);
    void fill(bool x, unsigned ix0, unsigned ix1);
    
    class Reader {
    protected:
      BitArray<W> *a;
      W w;
      unsigned ix;
    public:
      Reader(const BitArray<W> *host, unsigned ix0=0);
      unsigned index() const { return ix; }
      template<unsigned n>
      W read(); // reading off the end is safe and will just return 0 bits
      void seek(unsigned ix);
    };
    
    class Writer: public Reader {
    private:
      Reader::seek;
    public:
      Writer(BitArray<W> *host, unsigned ix0=0);
      ~Writer();
      template<unsigned n>
      void write(W x); // writing off the end is undefined!
      void flush();
    };
  };

  template<class W>
  class FastBitArray {
    static const unsigned logc = 9, bitc = 1u<<logc;
    Ref<BitArray<W> > bitsref;
    BitArray<W> *bits;
    Ref<unsigned> chksref;
    unsigned *chks;
  public:
    FastBitArray(unsigned len);
    class Builder {
      Ref<FastBitArray<W> > ref;
      typename BitArray<W>::Writer w;
      unsigned *pchk, chkpop;
    public:
      Builder(unsigned len);
      unsigned index() const { return w.index(); }
      template<unsigned n>
      void write(W x);
      Ref<FastBitArray<W> > finish();
    };
  public:
    const BitArray<W>* bit_array() const { return bits; }
    unsigned length() const { return bits->length(); }
    bool get(unsigned ix) const { return bits->get(ix); }
    unsigned count(unsigned ix0, unsigned ix1) const;
    unsigned find(unsigned ix0, unsigned nth) const;
  };
}
#endif
