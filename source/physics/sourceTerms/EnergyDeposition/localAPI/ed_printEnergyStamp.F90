!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_printEnergyStamp
!!
!! NAME
!!
!!  ed_printEnergyStamp
!!
!! SYNOPSIS
!!
!!  call ed_printEnergyStamp (real, intent (in) :: timeStep,
!!                            real, intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Prints out the current energy statistics of the laser:
!!
!!        i) total laser energy pumped into the domain so far
!!       ii) total laser energy exiting the domain unused so far
!!      iii) laser energy pumped into the domain at the current timestep
!!       iv) laser energy exiting the domain unused at the current timestep
!!
!! ARGUMENTS
!!
!!  timeStep       : Current time step value
!!  timeSimulation : Current simulation time
!!
!! NOTES
!!
!!  Only the master processor writes to the energy statistics file.
!!
!!***

subroutine ed_printEnergyStamp (timeStep, timeSimulation)

  implicit none
   
  real, intent (in) :: timeStep
  real, intent (in) :: timeSimulation

  return
end subroutine ed_printEnergyStamp
