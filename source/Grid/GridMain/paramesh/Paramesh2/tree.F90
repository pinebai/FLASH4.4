!!****ih* source/Grid/GridMain/paramesh/Paramesh2/tree
!!
!! NAME
!!
!!    tree
!!
!!
!! SYNOPSIS
!!
!!   data structure for storing tree information
!!   
!!
!! DESCRIPTION
!!
!!
!!
!! This is the include file for a quad or oct-tree data structure
!! The tree organizes a set of up to maxblocks_tr grids on each processor.
!! All the grids are assumed to be cartesian with a uniform size. Each 
!! grid has a level of refinement associated with it. The set of level 0
!! grids cover the computational domain without overlap. Each grid
!! can be the parent of 2**d offspring grids which completely cover 
!! the sub-domain of their parents, where d is the physical dimension
!! of the simulation. The linear resolution varies by a factor of 2 
!! between successive levels of refinement. At no point do we allow the
!! resolution to jump by more than one level of refinement.
!!
!!
!! In the following list the index i ranges from 1 to maxblocks. 
!!
!!       neigh(2,nfaces,i)     local and processor ids of block i's neighbors,
!!                               at i's refinement level. If a neighbor does 
!!                               not exist both values are set to -1, unless 
!!                               that face is at an external domain boundary
!!                               where non-periodic boundary conditions are to
!!                               be applied, in which case these are set to -20
!!                               or less, depending on the boundary conditions
!!                               to be applied on the boundary in question.
!!       child(2,nchild,i)     local and processor ids of block i's children
!!       parent(2,i)           local and processor ids of block i's parent
!!       coord(ndim,i)         array storing x,y and z coordinates of the
!!                               center of block i.
!!       bnd_box(2,ndim,i)     bounding box information for block i. The 
!!                               lower edge of block i along the j-th coordinate
!!                               axis is at bnd_box(1,j,i) and the upper edge
!!                               at bnd_box(2,j,i).
!!       bsize(ndim,i)          size of block i in the x, y and z directions.
!!       lrefine(i)            refinement level of block i.
!!       nodetype(i)           defines the node type, if 1 then the node is a
!!                               leaf node, if 2 then the node is a parent but
!!                               with at least 1 leaf child, otherwise it is
!!                               set to 3 and it does not have any up-to-date
!!                               data.
!!       empty(i)              used to designate empty blocks, for example
!!                               when an obstacle is inserted inside the
!!                               computational domain. normal blocks have
!!                               empty=0, empty blocks have empty=1.
!!       
!!       new_child(i)          if true then child has just been produced by
!!                               a refinement step, otherwise false.
!!       lnblocks              number of blocks on the local processor
!!       new_lnblocks          the new number of blocks on the local 
!!                               processor after a refinement or derefinement 
!!                               step.
!!       refine(i)             refinement flag. If set to .true. block i
!!                               will be refined during the next call to
!!                               REFINE_DEREFINE.
!!       derefine(i)           derefinement flag. If set to .true. block i
!!                               will be derefined during the next call to
!!                               REFINE_DEREFINE, provided this blocks parent
!!                               is not marked for refinement.
!! ADDED FOR EASIER AND MORE EFFICIENT MPI MESSAGING (KMO)
!!       neigh_type(nfaces,i)  types of the neighbors of block i
!!       child_type(nchild,i)  types of the children of block i
!!
!!
!!-----------------------------------------------------------------
!!
!!
!! ARGUMENTS
!!
!!***




module tree
  use physicaldata  
!
!
!
! $RCSfile: tree.F90,v $
! $Revision: 1.2 $
! $Date: 2004/09/07 00:58:09 $
!
!
!
  integer,save :: lrefine_min, lrefine_max
  integer maxblocks_tr
  parameter(maxblocks_tr=10*maxblocks)
!
! Number of children of a node
  integer nchild
  parameter(nchild=2**ndim)
!
! Number of faces on a grid block
  integer nfaces
  parameter(nfaces=2*ndim)
!
! Parameters used to define array sizes
  integer mdim,mchild,mfaces
  parameter(mdim=3,mchild=2**mdim,mfaces=2*mdim)
!
! Common block storing tree datastructure

!
  integer,target,dimension(2,mfaces,maxblocks_tr) :: neigh
  integer,target,dimension(2,mchild,maxblocks_tr) :: child
  integer,target,dimension(2,maxblocks_tr) :: parent
  integer,target,dimension(maxblocks_tr) :: lrefine

  integer lnblocks,new_lnblocks
  integer :: grid_changed = 1 ! added to support BH tree gravity implementation
!
  integer,target,dimension(maxblocks_tr) :: nodetype
  integer,target,dimension(mfaces,maxblocks_tr) :: neigh_type
  integer,target,dimension(mchild,maxblocks_tr) :: child_type
  integer,target,dimension(maxblocks_tr) :: empty
  
  logical,target,dimension(maxblocks_tr) :: newchild
  logical,target,dimension(maxblocks_tr) :: derefine
  logical,target,dimension(maxblocks_tr) :: refine
  logical,target,dimension(maxblocks_tr) :: stay

!
  real,target,dimension(maxblocks_tr) :: work_block
  real,target,dimension(mdim,maxblocks_tr) :: coord
  real,target,dimension(mdim,maxblocks_tr) :: bsize
  real,target,dimension(2,mdim,maxblocks_tr) :: bnd_box
!
!--------------------------------------------
!
! A convention is established for numbering the neighbors (or faces
! of a block. The first neighbor is at lower x coordinate, the 
! second at higher x, the third at lower y, fourth at higher y, fifth
! at lower z and the sixth at higher z.
!
! The convention by which the children of a block are numbered is the
! same as the fortran array ordering, so that the first child is
! at lower x, y and z coordinate, the second child is at upper x
! but lower y and z, the third is at lower x, upper y and lower z,
! and so on.
!
! When a block has a refined neighbor we will need to know which children
! of this neighbor are to provide guard cell information. The id's of the
! correct children are stored in kchild using the conventions described 
! above. For example, if we are working on the 3rd neighbor of the
! current block and it is at finer refinement level, then we must access
! the children designated by kchild(:,3), in this case children 1, 2, 5
! and 6.
!
!--------------------------------------------

end module tree

