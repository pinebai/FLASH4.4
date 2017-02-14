!!****f* source/Simulation/Simulation_computeAnalytical
!!
!! NAME
!!
!!  Simulation_computeAnalytical
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_computeAnalytical(integer(IN) :: blockID, 
!!                                    real(IN)    :: tcurr)
!!
!!
!!
!! DESCRIPTION
!!
!!  Compute an analytical solution.
!!
!!  This is simulation-dependent, there is no general implementation.
!! 
!! ARGUMENTS
!!
!!  blockID -        the number of the block to initialize
!!  tcurr   -        current time
!!
!! SIDE EFFECTS
!!
!!  The analytical solution is computed and stored in an appropriate slot
!!  (or slots) in the solution vector, UNK.
!!
!! NOTES
!! 
!!  By default this is just a stub that doesn't do anything.
!!
!!***

subroutine Simulation_computeAnalytical(blockId, tcurr)

  implicit none

  integer, intent(in) :: blockId
  real,    intent(in) :: tcurr

end subroutine Simulation_computeAnalytical
