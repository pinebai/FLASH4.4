!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhBuildTreeBlock
!!
!! NAME
!!
!!  gr_bhBuildTree
!!
!!
!! SYNOPSIS
!!
!!   gr_bhBuildTreeBlock(
!!        integer,intent(in) :: block
!!        )
!!
!! DESCRIPTION
!!
!!   Build a tree for a specific block.
!!
!! ARGUMENTS
!!
!!  block : ID of a block where the block-tree is constructed
!!
!!***


subroutine gr_bhBuildTreeBlock(block)

  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getCellCoords, &
    Grid_getBlkPtr, Grid_getDeltas, Grid_releaseBlkPtr
  use gr_bhLocalInterface, ONLY : gr_bhGetTreeSize, gr_bhFillBotNode &
  &   , gr_bhAccBotNode, gr_bhNormalizeNode, gr_bhAccNode
  use gr_bhLocalInterface, ONLY : gr_bhGetTreePos
  use gr_bhData, ONLY: myPE => gr_bhTreeMyPE, gr_bhTreeLevels, gr_bhTreeArray, &
       gr_bhLocCoords, gr_bhTreeBS, gr_bhTreeParentTree, GR_TREE_NSIZE, GR_TREE_BNSIZE, &
       GR_TREE_IX, GR_TREE_IY, GR_TREE_IZ, gr_bhTreeLrefine, gr_bhTreeBCen, gr_bhTreeCellSize
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"

  integer,intent(in) :: block
  integer       :: i, j, k, l, pos, fac, istat, level
  integer, dimension(2,MDIM)   :: blkLimits,blkLimitsGC
  real, DIMENSION(:,:,:,:), POINTER :: solnData
  integer, DIMENSION(MDIM) :: point
  integer       :: multi(1:gr_bhTreeLevels), mm
  real          :: del(MDIM)
  real          :: dvol, mass
  logical       :: gcell = .false.
  real          :: m = 0, xmc = 0, ymc = 0, zmc = 0
  real, allocatable :: botnode(:), accnode(:), node(:)
  

  ! allocate the tree array
  nullify(gr_bhTreeArray(myPE,block)%p)
  allocate(gr_bhTreeArray(myPE,block)%p(1:gr_bhGetTreeSize(gr_bhTreeLevels)), stat=istat)
  if (istat /= 0) then
    call Driver_abortFlash("could not allocate tree in gr_bhBuildTreeBlock.F90")
  endif

  ! fill tree with zeros
  do i = 1, gr_bhGetTreeSize(gr_bhTreeLevels)
    gr_bhTreeArray(myPE,block)%p(i) = 0
  enddo

  ! get information about coordinates of the block
  call Grid_getBlkIndexLimits(block,blkLimits,blkLimitsGC)
  call Grid_getCellCoords(IAXIS, block, CENTER, gcell, gr_bhLocCoords(:, IAXIS, block), gr_bhTreeBS)
  call Grid_getCellCoords(JAXIS, block, CENTER, gcell, gr_bhLocCoords(:, JAXIS, block), gr_bhTreeBS)
  call Grid_getCellCoords(KAXIS, block, CENTER, gcell, gr_bhLocCoords(:, KAXIS, block), gr_bhTreeBS)
  
  ! create the lowest level (and prepare the second lowest)
  call Grid_getBlkPtr(block,solnData,CENTER)
  call Grid_getDeltas(block,del)
  allocate(botnode(GR_TREE_BNSIZE), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate botnode in gr_bhBuildTreeBlock.F90")
  allocate(accnode(GR_TREE_NSIZE),  stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate accnode in gr_bhBuildTreeBlock.F90")
  do k = 1,gr_bhTreeBS
    point(KAXIS) = k+blkLimits(LOW,KAXIS)-1
    do j = 1,gr_bhTreeBS
      point(JAXIS) = j+blkLimits(LOW,JAXIS)-1
      do i = 1,gr_bhTreeBS
        point(IAXIS) = i+blkLimits(LOW,IAXIS)-1

        ! determine the multiindex
        do l = 1, gr_bhTreeLevels
          fac = 2**(gr_bhTreeLevels - l)
          multi(l) = 1 + mod((i-1),2*fac)/fac  &
                   + 2 *(mod((j-1),2*fac)/fac) &
                   + 4 *(mod((k-1),2*fac)/fac)
        enddo

        ! the lowest level - just write masses into the tree
        pos = gr_bhGetTreePos(gr_bhTreeLevels, multi)
        call gr_bhFillBotNode(block, point, blkLimits, solnData, botnode)
        gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_BNSIZE-1) = botnode

        ! prepare the second lowest level (sum mass and m*(x,y,z) to the tree array)
        multi(gr_bhTreeLevels) = 0
        pos = gr_bhGetTreePos(gr_bhTreeLevels-1, multi)
        accnode = gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_NSIZE-1)
        call gr_bhAccBotNode(block, point, blkLimits, gr_bhLocCoords, solnData, botnode, accnode)
        gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_NSIZE-1) = accnode
      enddo
    enddo
  enddo
  deallocate(botnode)

  ! second lowest level - normalization of the mass centers by the total mass is necessary
  allocate(node(GR_TREE_NSIZE),  stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate node in gr_bhBuildTreeBlock.F90")

  do l = 1,gr_bhTreeLevels-1
    multi(l) = 1
  enddo
  multi(gr_bhTreeLevels) = 0
  do
    ! normalize mass center by the total mass in the branch
    pos = gr_bhGetTreePos(gr_bhTreeLevels-1, multi)
    node = gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_NSIZE-1)
    call gr_bhNormalizeNode(node)
    gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_NSIZE-1) = node
    if (abs(node(GR_TREE_IX) - gr_bhTreeBCen(IAXIS, block, myPE)) &
    &  > 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), IAXIS)*gr_bhTreeBS) then
      print *, "XMC out of node!", block, myPE, node(GR_TREE_IX), gr_bhTreeBCen(IAXIS, block, myPE) &
      & , gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), IAXIS)*gr_bhTreeBS
    endif
    if (abs(node(GR_TREE_IY) - gr_bhTreeBCen(JAXIS, block, myPE)) &
    &  > 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), JAXIS)*gr_bhTreeBS) then
      print *, "YMC out of node!", block, myPE, node(GR_TREE_IY), gr_bhTreeBCen(JAXIS, block, myPE) &
      & , gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), JAXIS)*gr_bhTreeBS
    endif
    if (abs(node(GR_TREE_IZ) - gr_bhTreeBCen(KAXIS, block, myPE)) &
    &  > 0.5*gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), KAXIS)*gr_bhTreeBS) then
      print *, "ZMC out of node!", block, myPE, node(GR_TREE_IZ), gr_bhTreeBCen(KAXIS, block, myPE) &
      & , gr_bhTreeCellSize(gr_bhTreeLrefine(block, myPE), KAXIS)*gr_bhTreeBS
    endif

    
    multi(gr_bhTreeLevels-1) = multi(gr_bhTreeLevels-1) + 1 ! update multi-index counter
    ! rotate multi-index counter
    do l = gr_bhTreeLevels-1,2,-1
      if (multi(l) > 8) then
        multi(l) = 1
        if (l > 1) multi(l-1) = multi(l-1) + 1
      endif
    enddo
    if (multi(1) > 8) exit 
  enddo
  

  ! create higher levels
  do level = gr_bhTreeLevels-1,1,-1
    
    ! set the initial multi-index
    do l = 1,gr_bhTreeLevels
      if (l <= level) then
        multi(l) = 1
      else
        multi(l) = 0
      endif
    enddo
    
    ! set accumulators to zero
    accnode(1:GR_TREE_NSIZE) = 0.0
    do
      ! sum mass and mass centres from the octet of branches 
      pos = gr_bhGetTreePos(level, multi)
      node = gr_bhTreeArray(myPE, block)%p(pos:pos+GR_TREE_NSIZE-1)
      call gr_bhAccNode(node, accnode)

      ! update multi-index counter
      multi(level) = multi(level) + 1

      if (multi(level) > 8) then
        ! write sum to the parent node
        mm = multi(level) ! backup multi(level)
        multi(level) = 0
        pos = gr_bhGetTreePos(level-1, multi)
        multi(level) = mm
        call gr_bhNormalizeNode(accnode)
        gr_bhTreeArray(myPE,block)%p(pos:pos+GR_TREE_NSIZE-1) = accnode

        ! set accumulators to zero
        accnode(1:GR_TREE_NSIZE) = 0.0

      endif

      ! rotate the multi-index counter
      do l = level,2,-1
        if (multi(l) > 8) then
          multi(l) = 1
          if (l > 1) multi(l-1) = multi(l-1) + 1
        endif
      enddo
      if (multi(1) > 8) exit 
     
    enddo
  enddo

  gr_bhTreeParentTree(1:GR_TREE_NSIZE, block,myPE) = gr_bhTreeArray(myPE,block)%p(1:GR_TREE_NSIZE)

  deallocate(accnode)
  deallocate(node)
  call Grid_releaseBlkPtr(block,solnData)

  return
end subroutine gr_bhBuildTreeBlock

