!!****if* source/Grid/localAPI/gr_hypreExchangeFacB
!!
!!  NAME 
!!
!!  gr_hypreExchangeFacB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreExchangeFacB (integer,intent(IN) :: iFactorB,
!!                             integer,intent(IN) :: blockCount,
!!                             integer,intent(IN) :: blockList (blockCount))
!!
!!
!!  DESCRIPTION 
!!   This routine helps exchange iFactorB at fine-coarse boundaries. If called
!!   in UG mode this rountine returns without any action. In AMR mode it uses (tricks) 
!!   Grid_conserveFluxes to succesfully communicate the values across processors. Doing
!!   this saves us from performing more complicated point-point communication.
!!
!!
!! ARGUMENTS
!!   iFactorB   : Factor from the diffusion equation.
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!  
!! NOTES
!!   Stub implementation.
!!   
!!
!!***

!!REORDER(4): solnVec


subroutine gr_hypreExchangeFacB (iFactorB, blockCount, blockList) 
  
  implicit none
  
  integer,intent(IN) :: iFactorB
  integer,intent(IN) :: blockCount
  integer,intent(IN) :: blockList (blockCount)
  
  
  return
  
end subroutine gr_hypreExchangeFacB
