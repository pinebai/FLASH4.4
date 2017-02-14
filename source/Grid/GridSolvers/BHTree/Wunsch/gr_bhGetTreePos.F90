!!****if* source/Grid/GridSolvers/BHTree/Wunsch/gr_bhGetTreePos
!!
!! NAME
!!
!!  gr_bhGetTreePos
!!
!!
!! SYNOPSIS
!!
!!   call gr_bhGetTreePos(
!!        integer(in) :: level,
!!        integer(in) :: mi(1:gr_bhTreeLevels)
!!        )
!!
!! DESCRIPTION
!!
!!   Returns position of the node with multi-index mi in the tree array.
!!
!! ARGUMENTS
!!
!!  level  - level of the tree node
!!  mi     - multi-index (indexes 1-8 at individual levels, or 0 if it is 
!!            a higher level node)
!!
!! RESULT
!!
!!   Position of the node in the tree array.
!!
!!***



integer function gr_bhGetTreePos(level, mi)
  use gr_bhData, ONLY : gr_bhTreeLevels, GR_TREE_BNSIZE, GR_TREE_NSIZE
  implicit none
  integer,intent(in) :: level, mi(1:gr_bhTreeLevels)
  integer            :: pos, fs, l

  ! field size: 1 for the lowest level (only masses), 
  !             4 for the higher levels (mass, position of mass centre)
  if (level == gr_bhTreeLevels) then
    fs = GR_TREE_BNSIZE
  else
    fs = GR_TREE_NSIZE
  endif
  
  pos = 1 + GR_TREE_NSIZE*(8**level - 1)/7 ! offset of the level
  do l = level,1,-1
    pos = pos + (mi(l)-1)*fs
    fs = fs * 8
  enddo
  gr_bhGetTreePos = pos
end function gr_bhGetTreePos


