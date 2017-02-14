!!****f* source/Driver/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN):: blockCount,
!!                     integer(IN):: blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash, we call Driver_sourceTerms, which then
!!  makes the calls to Cool, Burn, Heat, and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the 
!!                 source term operator(s)
!!  dt           : the current timestep
!!
!!***

subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)

  implicit none

  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList
  integer, OPTIONAL, intent(IN):: pass

  return
end subroutine Driver_sourceTerms
