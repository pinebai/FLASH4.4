!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_traceRays
!!
!! NAME
!!
!!  ed_traceRays
!!
!! SYNOPSIS
!!
!!  call ed_traceRays (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Follows the collection of all rays at this stage through one block and deposits their energy
!!  in those blocks. This routine calls the appropriate subroutines according to the domain
!!  grid geometry specified.
!!
!! ARGUMENTS
!!
!!  timeStep : Current time step
!!
!!***

subroutine ed_traceRays (timeStep)

  implicit none

  real, intent (in) :: timeStep

  return
end subroutine ed_traceRays
