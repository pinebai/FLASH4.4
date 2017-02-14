#include "bittree_bitarray.hxx"
#include "bittree_mortontree.hxx"
#include "bittree_ref.hxx"

#include "mpi.h"

using namespace BitTree;

typedef unsigned W;

namespace { // private globals
  struct TheTree_ {
    virtual unsigned block_count() = 0;
    virtual void morton(int *lev, int *ijk, int *mort) = 0;
    virtual void refine_init() = 0;
    virtual void refine_mark(int lev, const int *ijk) = 0;
    virtual void refine_apply(MPI_Comm comm) = 0;
  };
  
  template<unsigned ndim>
  class TheTree: public TheTree_ {
    Ref<MortonTree<ndim,W> > tree;
    Ref<BitArray<W> > refine_delta;
  public:
    TheTree(const unsigned top[], const bool includes[]);
    unsigned block_count();
    void morton(int *lev, int *ijk, int *mort);
    void refine_init();
    void refine_mark(int lev, const int *ijk);
    void refine_apply(MPI_Comm comm);
  };
  
  Ref<TheTree_> the_tree;
}

template<unsigned ndim>
TheTree<ndim>::TheTree(const unsigned top[], const bool includes[]):
  tree(MortonTree<ndim,W>::make(top, includes)) {
}

template<unsigned ndim>
unsigned TheTree<ndim>::block_count() {
  return tree->blocks();
}

template<unsigned ndim>
void TheTree<ndim>::morton(int *lev, int *ijk, int *mort) {
  unsigned coord[ndim];
  for(unsigned d=0; d < ndim; d++)
    coord[d] = ijk[d];
  
  if(tree->inside(*lev, coord)) {
    typename MortonTree<ndim,W>::template Block<unsigned> b = tree->identify(*lev, coord);
    *lev = b.level;
    ijk[0] = b.coord[0];
    *mort = b.mort;
  }
  else {
    *lev = -1;
    *mort = -1;
  }
}

template<unsigned ndim>
void TheTree<ndim>::refine_init() {
  unsigned nbits = tree->id_upper_bound();
  refine_delta = BitArray<W>::make(nbits);
  refine_delta->fill(false);
}

template<unsigned ndim>
void TheTree<ndim>::refine_mark(
    int lev, // in: 0-based
    const int *ijk // in: 0-based
  ) {
  unsigned coord[ndim];
  for(unsigned d=0; d < ndim; d++)
    coord[d] = ijk[d];
  
  unsigned id = tree->identify(lev, coord).id;
  refine_delta->set(id, true);
}

template<unsigned ndim>
void TheTree<ndim>::refine_apply(MPI_Comm comm) {
  MPI_Allreduce(
    MPI_IN_PLACE,
    refine_delta->word_buf(),
    refine_delta->word_count(),
    MPI_UNSIGNED,
    MPI_BOR,
    comm
  );
  tree = tree->refine(refine_delta);
  refine_delta.nullify();
}


extern "C" bool bittree_initialized() {
  return !!the_tree;
}

extern "C" void bittree_init(
    const int *ndim,
    const int topsize[], // in: 0-based
    const bool includes[] // in: includes[topsize[ndim-1]]...[topsize[0]]
  ) {
  unsigned top[3];
  for(unsigned d=0; d < *ndim; d++)
    top[d] = topsize[d];
  
  switch(*ndim) {
  case 1: {
    Ref<TheTree<1> > r; new(r.alloc()) TheTree<1>(top, includes);
    the_tree = r;
    } break;
  case 2: {
    Ref<TheTree<2> > r; new(r.alloc()) TheTree<2>(top, includes);
    the_tree = r;
    } break;
  case 3: {
    Ref<TheTree<3> > r; new(r.alloc()) TheTree<3>(top, includes);
    the_tree = r;
    } break;
  }
}

extern "C" void bittree_block_count(int *count) {
  *count = the_tree->block_count();
}

// given a level and ijk block position for that level, find the finest block
// that exactly matches or just contains that block position
extern "C" void bittree_morton(
    int *lev, // inout: 0-based
    int *ijk, // inout: 0-based
    int *mort // out: 0-based
  ) {
  the_tree->morton(lev, ijk, mort);
}

// how to refine/derefine a tree:
//   call bittree_refine_init
//   call bittree_refine_mark for every local block thats going to change
//   call bittree_refine_apply
// after that the tree corresponds to the new domain
extern "C" void bittree_refine_init() {
  if(!!the_tree)
    the_tree->refine_init();
}

// mark block for nodetype change
extern "C" void bittree_refine_mark(
    const int *lev, // in: 0-based
    const int *ijk // in: 0-based
  ) {
  if(!!the_tree)
    the_tree->refine_mark(*lev, ijk);
}

extern "C" void bittree_refine_apply(int *comm_) {
  if(!!the_tree) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_);
    the_tree->refine_apply(comm);
  }
}

#if 0
int main() {
  unsigned size[2] = {3,3};
  unsigned excl[1][2] = {1,1};
  Ref<MortonTree<2,unsigned> > tree = MortonTree<2,unsigned>::make(size, 1, excl);
  //for(int i=0; i < 10; i++) {
  for(int i=0; i < 3; i++) {
    Ref<BitArray<unsigned> > delta = BitArray<unsigned>::make(tree->level_id0(0) + tree->blocks());
    delta->fill(false, 0, tree->level_id0(tree->levels()-1));
    delta->fill(true, tree->level_id0(tree->levels()-1), tree->id_upper_bound());
    //delta->set(tree->level_id0(tree->levels()-1), true);
    //delta->set(tree->blocks()-1, true);
    tree = tree->refine(delta);
    cout << "i=" << i << " tree blocks=" << tree->blocks() << '\n';
  }
  if(false) {
    unsigned sum = 0;
    for(int i=0; i < 1<<14; i++) {
      unsigned x[2] = { i%(1<<11), i%(1<<11) };
      MortonTree<2,unsigned>::Block<unsigned> b0 = tree->identify(10, x);
      /*MortonTree<2,unsigned>::Block<unsigned> b1 = tree->locate<unsigned>(b0.id);
      DBG_ASSERT(b0.id == b1.id);
      DBG_ASSERT(b0.level == b1.level);
      DBG_ASSERT(b0.mort_ix == b1.mort_ix);
      DBG_ASSERT(b0.coord[0] == b1.coord[0] && b0.coord[1] == b1.coord[1]);
      */sum = 5*sum ^ b0.mort;
    }
    cout << "sum:" << sum << '\n';
  }
  else
    cout << (MortonTree<2,unsigned>*)tree;
}
#endif
