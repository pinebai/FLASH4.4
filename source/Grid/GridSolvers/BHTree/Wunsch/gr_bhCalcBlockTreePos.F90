!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhCalcBlockTreePos
!!
!! NAME
!!
!!  gr_bhCalcBlockTreePos
!!
!!
!! SYNOPSIS
!!
!!  call gr_bhCalcBlockTreePos(
!!        )
!!
!! DESCRIPTION
!!
!!   Calculates positions of all nodes in the tree array and stores them in
!!   the gr_bhBlockTreePos array. Calculates also geometric centres of nodes
!!   in that array (with respect to the block centre, meassured in grid cells).
!!
!! ARGUMENTS
!!
!!
!! RESULT
!!
!!   Calculates the gr_bhBlockTreePos array:
!!     1. index: 
!!        GR_TREE_BTP_POS:  position in block tree
!!        GR_TREE_BTP_LEV:  level in block tree
!!        GR_TREE_BTP_IJK:  if bottom level: cell indeces, otherwise -1, -1, -1
!!        GR_TREE_BTP_C1-8: position of children in gr_bhBlockTreePos; if bottom level: -1 ... -1
!!
!!***



subroutine gr_bhCalcBlockTreePos()
  use gr_bhData, ONLY : gr_bhTreeLevels, gr_bhBlockTreePos, gr_bhTreeBS, &
      gr_bhTreeLoff, GR_TREE_BNSIZE, GR_TREE_NSIZE, GR_TREE_BTP_N, GR_TREE_BTP_POS, &
      GR_TREE_BTP_LEV, GR_TREE_BTP_I, GR_TREE_BTP_J, GR_TREE_BTP_K, &
      GR_TREE_BTP_C1, GR_TREE_BTP_C8, gr_bhBlockTreeNodeCen
  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"
  integer :: i, j, k, l, sp, level, istat
  integer :: multi(1:gr_bhTreeLevels), mi(1:gr_bhTreeLevels)
  integer :: fs, fac, pos, btp, cbtp
  integer :: stack(1:gr_bhTreeLevels, 1:(8**(gr_bhTreeLevels+1)-1)/7)

  ! allocate BlockTreePos
  allocate(gr_bhBlockTreePos(GR_TREE_BTP_N, (8**(gr_bhTreeLevels+1) - 1)/7), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhBlockTreePos in gr_bhCalcBlockTreePos")
  gr_bhBlockTreePos(:,:) = -1
  ! allocate BlockTreeNodeCen
  allocate(gr_bhBlockTreeNodeCen(MDIM, (8**(gr_bhTreeLevels+1) - 1)/7), stat=istat)
  if (istat /= 0) call Driver_abortFlash("could not allocate gr_bhBlockTreeNodeCen in gr_bhCalcBlockTreePos")

  ! put root node on stack
  do l = 1, gr_bhTreeLevels
    multi(l) = 0
  enddo
  stack(:,1) = multi(:)
  sp = 1

  ! Set the centre of the block into the BlockTreeNodeCen array
  gr_bhBlockTreeNodeCen(:,:) = -1000.0
  gr_bhBlockTreeNodeCen(:,1) = 0.0
  
  ! walk through the tree
  do
    ! take the node on the bottom of the stack
    multi(:) = stack(:,sp)
    sp = sp - 1

    ! determine level of the multi-index
    level = gr_bhTreeLevels
    do l = 1,gr_bhTreeLevels
      if (multi(l) == 0) then
        level = l - 1
        exit
      endif
    enddo

    ! compute position in the gr_bhBlockTreePos array
    fs = 1
    btp = 1 + (8**level - 1)/7 ! offset of the level
    do l = level,1,-1
      btp = btp + (multi(l)-1)*fs
      fs = fs * 8
    enddo

    ! compute position in the tree_array
    if (level == gr_bhTreeLevels) then
      fs = GR_TREE_BNSIZE
    else
      fs = GR_TREE_NSIZE
    endif
    pos = gr_bhTreeLoff(level)
    do l = level,1,-1
      pos = pos + (multi(l)-1)*fs
      fs = fs * 8
    enddo
    gr_bhBlockTreePos(GR_TREE_BTP_POS, btp) = pos
    gr_bhBlockTreePos(GR_TREE_BTP_LEV, btp) = level

    ! convert multi-index to indeces in the block
    if (level == gr_bhTreeLevels) then
      i = 1
      j = 1
      k = 1
      do l = 1,gr_bhTreeLevels
        fac = 2**(gr_bhTreeLevels-l)
        i = i + (mod(multi(l)-1,2))      * fac
        j = j + (mod((multi(l)-1)/2,2))  * fac
        k = k + (mod((multi(l)-1)/4,2))  * fac
      enddo
      gr_bhBlockTreePos(GR_TREE_BTP_I, btp) = i
      gr_bhBlockTreePos(GR_TREE_BTP_J, btp) = j
      gr_bhBlockTreePos(GR_TREE_BTP_K, btp) = k
      !print *, "BTP NodeCen: ", i,j,k, gr_bhBlockTreeNodeCen(:, btp)
    else

      ! put all children on the stack
      mi(:) = multi(:)
      do i = 1,8
        mi(level+1) = i
        sp = sp + 1
        stack(:,sp) = mi(:)
        ! determine position of the child in the BlockTreePos array
        fs = 1
        cbtp = 1 + (8**(level+1) - 1)/7 ! offset of the level
        do l = level+1,1,-1
          cbtp = cbtp + (mi(l)-1)*fs
          fs = fs * 8
        enddo
        gr_bhBlockTreePos(GR_TREE_BTP_C1-1+i, btp) = cbtp

        ! set node centres of children
        gr_bhBlockTreeNodeCen(IAXIS, cbtp) = gr_bhBlockTreeNodeCen(IAXIS, btp) &
        & + (gr_bhTreeBS/2**(level+1))*(mod(i-1,2) - 0.5)
        gr_bhBlockTreeNodeCen(JAXIS, cbtp) = gr_bhBlockTreeNodeCen(JAXIS, btp) &
        & + (gr_bhTreeBS/2**(level+1))*(mod((i-1)/2,2) - 0.5)
        gr_bhBlockTreeNodeCen(KAXIS, cbtp) = gr_bhBlockTreeNodeCen(KAXIS, btp) &
        & + (gr_bhTreeBS/2**(level+1))*(mod((i-1)/4,2) - 0.5)

      enddo
    endif
              
    ! stack is empty - exiting
    if (sp == 0) exit
  enddo

  !do i = 1, (8**(gr_bhTreeLevels+1)-1)/7
  !    print *, "BTP: ", i, gr_bhBlockTreePos(:, i)
  !enddo
  return
end subroutine gr_bhCalcBlockTreePos



