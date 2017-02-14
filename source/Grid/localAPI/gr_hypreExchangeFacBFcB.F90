!!****if* source/Grid/localAPI/gr_hypreExchangeFacBFcB
!!
!!  NAME 
!!
!!  gr_hypreExchangeFacBFcB
!!
!!  SYNOPSIS
!!
!!  call gr_hypreExchangeFacBFcB (integer,intent(IN) :: iFactorB,
!!                             integer,intent(IN) :: blockCount,
!!                             integer, dimension(blockCount),intent(IN):: blockList)
!!
!!  DESCRIPTION 
!!   This routine helps exchange iFactorB at fine-coarse boundaries. If called
!!   in UG mode this rountine returns without any action. In AMR mode it uses (tricks) 
!!   Grid_conserveFluxes to succesfully communicate the values across processors. Doing
!!   this saves us from performing more complicated point-point communication.
!!
!!   This FcB ("face-centered B") variant expects iFactorB to be an index to
!!   a variable in an allocated face-centered scratch (from Grid_ascModule).
!!
!! ARGUMENTS
!!   iFactorB   : Factor from the diffusion equation.
!!   blockCount : The number of blocks in the list.   
!!   blockList  : The list of blocks on which the solution must be updated.  
!!
!!
!! SIDE EFFECTS
!!
!!   Upon return, data will have been stored (and corrected) in the global storage
!!   area that he Grid unit uses for storage of flux data.  The data can be retrieved
!!   by calling Grid_getFluxData.
!!
!! NOTES
!!
!!   With UG, there is only this stub.
!!
!! SEE ALSO
!!
!!   Grid_ascModule
!!
!!***

subroutine gr_hypreExchangeFacBFcB (iFactorB, blockCount, blockList)
  implicit none
  integer,intent(IN) :: iFactorB
  integer,intent(IN) :: blockCount
  integer, dimension(blockCount),intent(IN):: blockList
end subroutine gr_hypreExchangeFacBFcB
