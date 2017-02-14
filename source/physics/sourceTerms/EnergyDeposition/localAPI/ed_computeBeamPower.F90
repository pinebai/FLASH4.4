!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_computeBeamPower
!!
!! NAME
!!
!!  ed_computeBeamPower
!!
!! SYNOPSIS
!!
!!  call ed_computeBeamPower (real,    intent (in)  :: timeStep,
!!                            real,    intent (in)  :: timeSimulation,
!!                            integer, intent (in)  :: pulseID, 
!!                            real,    intent (out) :: beamPower)
!!
!! DESCRIPTION
!!
!!  Computes the average power of the beam between the time frame:
!!
!!             timeSimulation --> timeSimulation + timeStep
!!
!!  from the pulse data supplied. Zero beam power is returned, if the time
!!  frame is outside the pulse defining time range.
!!
!! ARGUMENTS
!!
!!  timeStep       : current time step value
!!  timeSimulation : current simulation time
!!  pulseID        : pulse identification number
!!  beamPower      : the average beam power obtained
!!
!! NOTES
!!
!!***

subroutine ed_computeBeamPower (timeStep, timeSimulation, pulseID, beamPower)

  implicit none

  real,    intent (in)  :: timeStep
  real,    intent (in)  :: timeSimulation
  integer, intent (in)  :: pulseID
  real,    intent (out) :: beamPower

  beamPower = 0.0

  return
end subroutine ed_computeBeamPower
