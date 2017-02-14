!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commIncrementDeadRays
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
  implicit none
  integer, intent(IN) :: numDeadRays
end subroutine ed_commIncrementDeadRays
