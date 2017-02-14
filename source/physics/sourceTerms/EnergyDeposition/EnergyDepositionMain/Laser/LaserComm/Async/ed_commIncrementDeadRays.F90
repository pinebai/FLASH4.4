!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commIncrementDeadRays
!!
!!  NAME     
!!   ed_commIncrementDeadRays
!!
!!  SYNOPSIS
!!   ed_commIncrementDeadRays(integer, intent(IN) :: numDeadRays)
!!
!!  DESCRIPTION
!!    Increments the count of "dead" rays.  These are rays which have been
!!    absorbed or have left the domain.
!!
!!  ARGUMENTS
!!    numDeadRays : Number of dead rays
!!
!!***

subroutine ed_commIncrementDeadRays(numDeadRays)
  use UTCounter_sharedCounter, ONLY : UTCounter_incrementCounter
  use Driver_interface, ONLY : Driver_abortFlash
  implicit none
  integer, intent(IN) :: numDeadRays
  if (numDeadRays < 0) then
     call Driver_abortFlash("Negative number of rays???")
  end if
  !$omp critical (IncrementDeadRays)
  call UTCounter_incrementCounter(numDeadRays)
  !$omp end critical (IncrementDeadRays)
end subroutine ed_commIncrementDeadRays
