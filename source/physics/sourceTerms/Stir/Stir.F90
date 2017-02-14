!!****f* source/physics/sourceTerms/Stir/Stir
!!
!! NAME
!!
!!  Stir
!!
!! SYNOPSIS
!!
!!  Stir(integer(IN)::blockCount
!!       integer(IN)::blockList(blockCount),
!!       real(IN) :: dt,
!!
!! DESCRIPTION
!! Apply the isothermal cooling and stirring operator 
!! on the list of blocks provided as input
!!
!! ARGUMENTS
!! blockCount   : The number of blocks in the list
!! blockList(:) : The list of blocks on which to apply the stirring operator
!! dt           : the current timestep
!!
!!***

subroutine Stir(blockCount,blockList,dt)
  implicit none
  integer,intent(IN) :: blockCount
  integer,dimension(blockCount), intent(IN) :: blockList
  real,intent(IN) :: dt

  return
end subroutine Stir
