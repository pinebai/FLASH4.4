!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!
!!****h* headers/tree
!!
!! NAME
!!
!!   physicaldata
!!
!! SYNOPSIS
!!
!!   module tree
!!
!! INCLUDES
!!
!!   paramesh_preprocessor.fh
!!
!! USES
!!
!!   paramesh_dimensions
!!
!! DESCRIPTION
!!
!!   Fortran 90 Module which 'holds' the tree data which PARAMESH uses
!!   for a quad or oct-tree data structure, implemented on a parallel computer.
!!   This data includes information such which blocks are neighbors of other 
!!   blocks, which block is another block's parent, and which blocks are child 
!!   blocks of a particular block.
!!
!!   A convention is established for numbering the neighbors (or faces
!!   of a block. The first neighbor is at lower x coordinate, the 
!!   second at higher x, the third at lower y, fourth at higher y, fifth
!!   at lower z and the sixth at higher z.
!!
!!   The convention by which the children of a block are numbered is the
!!   same as the fortran array ordering, so that the first child is
!!   at lower x, y and z coordinate, the second child is at upper x
!!   but lower y and z, the third is at lower x, upper y and lower z,
!!   and so on.
!!
!!   When a block has a refined neighbor we will need to know which children
!!   of this neighbor are to provide guard cell information. The id's of the
!!   correct children are stored in kchild using the conventions described 
!!   above. For example, if we are working on the 3rd neighbor of the
!!   current block and it is at finer refinement level, then we must access
!!   the children designated by kchild(:,3), in this case children 1, 2, 5
!!   and 6.
!!
!!   The tree organizes a set of up to maxblocks_tr grids on each processor.
!!   All the grids are assumed to be cartesian with a uniform size. Each 
!!   grid has a level of refinement associated with it. The set of level 0
!!   grids cover the computational domain without overlap. Each grid
!!   can be the parent of 2**d offspring grids which completely cover 
!!   the sub-domain of their parents, where d is the physical dimension
!!   of the simulation. The linear resolution varies by a factor of 2 
!!   between successive levels of refinement. At no point do we allow the
!!   resolution to jump by more than one level of refinement.
!!
!!
!!
!! AUTHORS
!!
!!  Peter MacNeice and Kevin Olson
!!
!!***

!!REORDER(5): unk, facevar[xyz], tfacevar[xyz]
!!REORDER(4): recvar[xyz]f
#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------
! tree module



      module tree

      use paramesh_dimensions

      private

      public :: maxblocks_tr
      public :: nchild, nfaces, mchild, mfaces, mdim,  mflags
      public :: mboundaries,nboundaries

! Block limit to be used in manipulating tree when modifying the grid.
#ifndef LIBRARY
      integer, parameter :: maxblocks_tr=10*maxblocks
#else
      integer :: maxblocks_tr
#endif
! Number of children of a node
#ifndef LIBRARY
      integer, parameter :: nchild=2**ndim
#else
      integer :: nchild
#endif

! Number of faces on a grid block
#ifndef LIBRARY
      integer, parameter :: nfaces=2*ndim
#else
      integer :: nfaces
#endif

! Parameters used to define array sizes
      integer, parameter :: mdim=3,mchild=2**mdim,mfaces=2*mdim

! Parameters used to declare the number of block marker flags needed
#ifndef LIBRARY
! --<< USER EDIT >>--
      integer, parameter :: mflags=MFLAGS
#else
      integer, save      :: mflags
#endif

! Parameters used to declare the number of boundary regions where boundary
! conditions are to be applied.  Typically: 2*ndim
#ifndef LIBRARY
! --<< USER EDIT >>--
      integer, parameter :: nboundaries=NBOUNDARIES
#else
      integer,save       :: nboundaries
#endif

      integer, parameter :: mboundaries=100

      public :: neigh,child,which_child
      public :: parent,lrefine,lnblocks,new_lnblocks
      public :: nodetype,empty,bflags,newchild,derefine,refine
      public :: stay,work_block,coord,bsize,bnd_box
      public :: grid_xmin,grid_xmax,grid_ymin,grid_ymax
      public :: grid_zmin,grid_zmax
      public :: lrefine_max,lrefine_min
      public :: level_cell_sizes

! Variables for storing tree datastructure
#ifndef LIBRARY
      integer, save :: neigh(2,mfaces,maxblocks_tr)
      integer, save :: child(2,mchild,maxblocks_tr)
      integer, target, save :: which_child(maxblocks_tr)
      integer, save :: parent(2,maxblocks_tr),lrefine(maxblocks_tr)
#else
      integer, allocatable, save :: neigh(:,:,:)
      integer, allocatable, save :: child(:,:,:)
      integer, target, allocatable, save :: which_child(:)
      integer, allocatable, save :: parent(:,:),lrefine(:)
#endif
      integer, save :: lnblocks,new_lnblocks
#ifndef LIBRARY
      integer, save :: nodetype(maxblocks_tr)
      integer, save :: empty(maxblocks_tr)
      integer, target, save :: bflags(mflags,maxblocks_tr)
      logical, save :: newchild(maxblocks_tr)
      logical, save :: derefine(maxblocks_tr),refine(maxblocks_tr)
      logical, save :: stay(maxblocks_tr)
#else
      integer, allocatable, save :: nodetype(:)
      integer, allocatable, save :: empty(:)
      integer, target, allocatable, save :: bflags(:,:)
      logical, allocatable, save :: newchild(:)
      logical, allocatable, save :: derefine(:),refine(:)
      logical, allocatable, save :: stay(:)
#endif
#ifndef LIBRARY
      real, save :: work_block(maxblocks_tr)
      real, save :: coord(mdim,maxblocks_tr)
      real, save :: bsize(mdim,maxblocks_tr)
      real, save :: bnd_box(2,mdim,maxblocks_tr)
#else
      real, allocatable, save :: work_block(:)
      real, allocatable, save :: coord(:,:)
      real, allocatable, save :: bsize(:,:)
      real, allocatable, save :: bnd_box(:,:,:)
#endif
      real,save :: grid_xmin,grid_xmax
      real,save :: grid_ymin,grid_ymax
      real,save :: grid_zmin,grid_zmax

#ifndef LIBRARY
      real,save :: level_cell_sizes(mdim,maxlevels)
#else
      real, allocatable, save :: level_cell_sizes(:,:)
#endif
      integer, save :: lrefine_max,lrefine_min

! flag to record grid change
      public :: grid_changed, grid_analysed_mpi
      integer, save :: grid_changed, grid_analysed_mpi

! added for surrblks calculation with mpi
      public :: boundary_box,boundary_index
#ifndef LIBRARY
      real, save    :: boundary_box(2,mdim,mboundaries)
      integer, save :: boundary_index(mboundaries)
#else
      real, allocatable,save    :: boundary_box(:,:,:)
      integer, allocatable, save :: boundary_index(:)
#endif

! added for use with mpi block buffering
      public :: strt_buffer,last_buffer
      public :: strt_buffer_tree,last_buffer_tree
      public :: laddress,surr_blks
#ifdef SAVE_MORTS
      public :: surr_morts
#endif
      integer, save :: strt_buffer,last_buffer
      integer, save :: strt_buffer_tree,last_buffer_tree
#ifndef LIBRARY
      integer, save :: surr_blks(3,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
#ifdef SAVE_MORTS
      integer, save :: surr_morts(6,3,1+2*k2d,1+2*k3d,maxblocks_alloc)
#endif
      integer, save :: laddress(1:2,1:maxblocks_alloc)
#else
      integer, allocatable, save :: surr_blks(:,:,:,:,:)
#ifdef SAVE_MORTS
      integer, allocatable, save :: surr_morts(:,:,:,:,:)
#endif
      integer, allocatable, save :: laddress(:,:)
#endif


! arrays to store info about block neighbors which are boundaries
      public :: bc_block_neighs,bc_block_neighs_send
      public :: bc_block_neighs_length
      public :: bc_block_neighs_status
      integer,save,allocatable :: bc_block_neighs(:,:)
      integer,save,allocatable :: bc_block_neighs_send(:,:)
      integer,save             :: bc_block_neighs_length
      integer,save             :: bc_block_neighs_status


! DECLARE variables which are targets
      target refine, derefine, newchild, empty
      target lrefine, nodetype, work_block
      target parent, coord, bsize, neigh
      target child, bnd_box, stay


!!****v* tree/neigh
!!
!! NAME
!!
!!   neigh
!!
!! SYNOPSIS
!!
!!   public, integer :: neigh(2,mfaces,maxblocks_tr)
!!
!!   public, integer, allocatable :: neigh(:,:,:)
!!
!! DESCRIPTION
!!
!!   Local and processor ids of a block's neighbors, at that block's same 
!!   refinement level.  If a neighbor does not exist both values are set to -1,
!!   unless that face is at an external domain boundary where non-periodic 
!!   boundary conditions are to be applied, in which case these are set to 
!!   -20 or less, depending on the boundary conditions to be applied on the 
!!   boundary in question.
!!
!!***

!!****v* tree/child
!!
!! NAME
!!
!!   child
!!
!! SYNOPSIS
!!
!!   public, integer :: child(2,mchild,maxblocks_tr)
!!
!!   public, integer, allocatable :: child(:,:,:)
!!
!! DESCRIPTION
!!
!!   Local and processor ids of a block's children.
!!
!!***

!!****v* tree/parent
!!
!! NAME
!! 
!!   parent
!!
!! SYNOPSIS
!!
!!   public, integer :: parent(2,maxblocks_tr)
!!
!!   public, integer, allocatable :: parent(:,:)
!!
!! DESCRIPTION
!!  
!!   Local and processor ids of a block's parent.
!!
!!***

!!****v* tree/coord
!!
!! NAME
!!
!!   coord
!!
!! SYNOSIS
!!
!!   public, real :: coord(mdim, maxblocks_tr)
!!
!!   public, real, allocatable :: coord(mdim, maxblocks_tr)
!!
!! DESCRIPTION
!!
!!   An array storing x,y and z coordinates of the center of a block.
!!
!!***

!!****v* tree/bnd_box
!!
!! NAME
!!
!!   bnd_box
!!
!! SYNOPSIS
!!
!!   public, real :: bnd_box(2,mdim,maxblocks_tr)
!!
!!   public, real, allocatable :: bnd_box(:,:,:)
!!
!! DESCRIPTION
!!
!!   The bounding box information for a block. The lower edge of block i along
!!   along the j-th coordinate axis is at bnd_box(1,j,i) and the upper edge
!!   is at bnd_box(2,j,i).
!!
!!***

!!****v* tree/bsize
!! NAME
!!
!!   bsize
!!
!! SYNOPSIS
!!
!!   public, real :: bsize(mdim,maxblocks_tr)
!!
!!   public, real, allocatable :: bsize(:,:)
!!
!! DESCRIPTION
!!
!!   Physical size of a block in the x, y and z directions.
!!
!!***

!!****v* tree/lrefine
!!
!! NAME
!!
!!   lrefine
!!
!! SYNOPSIS
!!
!!   public, integer :: lrefine(maxblocks_tr)
!!
!!   public, integer, allocatable :: lrefine(:)
!!
!! DESCRIPTION
!!
!!   The refinement level of a block.
!!
!!***

!!****v* tree/nodetype
!!
!! NAME
!!
!!   nodetype
!!
!! SYNOPSIS
!!
!!   public, integer :: nodetype(maxblocks_tr)
!!   
!!   public, integer, allocatable :: nodetype(:)
!!
!! DESCRIPTION
!!
!!   Defines the node type, if 1 then the node is a leaf node, if 2 then the 
!!   node is a leaf node, if 2 then the node is a parent but with at least
!!   1 leaf child, otherwise it is set to 3 and it does not have any up-to-date
!!   data.
!!
!!***

!!****v* tree/empty
!!
!! NAME
!!
!!    empty
!!
!! SYNOPSIS
!!
!!   public, integer :: empty(maxblocks_tr)
!!
!!   public, integer, allocatable :: empty(:)
!!
!! DESCRIPTION
!!
!!   Used to designate empty blocks, for example, when an obstacle is inserted 
!!   inside the computational domain. normal blocks have empty=0, empty blocks 
!!   have empty=1.
!!
!!***

!!****v* tree/bflags
!!
!! NAME
!!
!!   bflags
!! 
!! SYNOPSIS
!!
!!   public, integer :: bflags(mflags,maxblocks_tr)
!!
!!   public, integer, allocatable :: bflags(:,:)
!!
!! DESCRIPTION
!!
!!   An array of integer flags which can be used to control computation on the 
!!   grid blocks and which are inherited by children from their parents.
!!
!!***

!!****v* tree/which_child
!!
!! NAME
!!
!!   which_child
!!
!! SYNOPSIS
!!
!!   public, integer :: which_child(maxblocks_tr)
!!
!!   public, integer, allocatable :: which_child(:)
!!
!! DESCRIPTION
!!
!!   An integer identifying which part of the parents volume this child 
!!   corresponds to.
!!
!!***

!!****v* tree/newchild
!!
!! NAME
!!
!! newchild
!!
!! SYNOPSIS
!!
!!   public, logical :: new_child(maxblocks_tr)
!!
!!   public, logical, allocatable :: newchild(:)
!!
!! DESCRIPTION
!!   
!!   If true then child has just been produced by a refinement step, otherwise 
!!   false.
!!
!!***

!!****v* tree/lnblocks
!!
!! NAME
!!
!!   lnblocks
!!
!! SYNOPSIS
!!
!!   public, integer :: lnblocks
!!
!! DESCRIPTION
!!
!!   The number of blocks on the local processor.
!!
!!***

!!****v* tree/new_lnblocks
!!
!! NAME
!!  
!!   new_lnblocks
!!
!! SYNOPSIS
!! 
!!   public, integer :: new_lnblocks
!!
!! DESCRIPTION
!!
!!   The new number of blocks on the local processor after a refinement or 
!!   derefinement step.
!!
!!***

!!****v* tree/refine
!!
!! NAME
!!
!!   refine
!!
!! SYNOPSIS
!!
!!   public, logical :: refine(maxblocks_tr)
!!
!!   public, logical, allocatable :: refine(:)
!!    
!! DESCRIPTION
!!
!!   The refinement flag. If set to .true. for a block, then that block will be
!!   refined during the next call to 'amr_refine_derefine'.
!!
!!***

!!****v* tree/derefine
!!
!! NAME
!!
!!   derefine
!!
!! SYNOPSIS
!!
!!   public, logical :: derefine(maxblocks_tr)
!!
!!   public, logical, allocatable :: derefine(:)
!!    
!! DESCRIPTION
!!
!!  The derefinement flag. If set to .true. for a block, then PARAMESH will 
!!  attempt to derefine that block at the next call to 'amr_refine_derefine' if 
!!  it comforms to the rules for derefinement (i.e. all siblings marked for 
!!  derefinement and will not produce a refinement jump of more than 1 level).
!!
!!***

!!****v* tree/stay
!!
!! NAME
!!
!!   stay
!!
!! SYNOPSIS
!!
!!   public, logical :: stay(maxblocks_tr)
!!
!!   public, logical, allocatable :: stay(:)
!!    
!! DESCRIPTION
!!
!!  If set to .true. for a block, then PARAMESH will neither refine or derefine
!!  the block at the next call to 'amr_refine_derefine'.
!!
!!***

!!****v* tree/work_block
!!
!! NAME
!!
!!   work_block
!!
!! SYNOPSIS
!!
!!   public, real :: work_block(maxblocks_tr)
!!
!!   public, real, allocatable :: work_block(:)
!!    
!! DESCRIPTION
!!
!!    The work weighting given to each block from which the load balancing 
!!    across processors is calculated.
!!
!!***

!-----------------------------------------------------------------
	
      end module tree
