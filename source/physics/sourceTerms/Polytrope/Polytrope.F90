!!****f* source/physics/sourceTerms/Polytrope/Polytrope
!!
!! NAME
!!  Polytrope
!!
!! SYNOPSIS
!!  Polytrope(integer(IN)::blockCount
!!            integer(IN)::blockList(blockCount),
!!            real(IN)::dt)
!!
!! DESCRIPTION
!!  Implement the polytropic eos as source term
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList(:) : The list of blocks on which to apply the Polytrope operator
!!  dt           : the current timestep
!!
!! WRITTEN BY
!!   Christoph Federrath 2007
!!
!!***

subroutine Polytrope(blockCount,blockList,dt)
  implicit none
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN) :: blockList
  real, intent(IN) :: dt
end subroutine Polytrope
