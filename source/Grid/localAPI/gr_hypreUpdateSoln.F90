!!****if* source/Grid/localAPI/gr_hypreUpdateSoln
!!
!!  NAME 
!!
!!  gr_hypreUpdateSoln
!!
!!  SYNOPSIS
!!
!!  call gr_hypreUpdateSoln (integer,intent(IN) :: iVar,
!!                           integer,intent(IN) :: blockCount,
!!                           integer,intent(IN) :: blockList (blockCount))
!!
!!  DESCRIPTION 
!!      This routine updates solution after solve (diffusion operation ,AX=B).
!!
!! ARGUMENTS
!!   iVar       : Variable on which the diffusion operatoion is performed (e.g TEMP_VAR)
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!
!!   Stub implementation.   
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreUpdateSoln (iVar, blockCount, blockList)
  
  implicit none

  integer,intent(IN) :: iVar
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList (blockCount)
  
  
  return
  
end subroutine gr_hypreUpdateSoln
