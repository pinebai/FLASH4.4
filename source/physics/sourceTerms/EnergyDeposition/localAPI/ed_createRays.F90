!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_createRays
!!
!! NAME
!!
!!  ed_createRays
!!
!! SYNOPSIS
!!
!!  call ed_createRays (integer, intent (in) :: blockCount, 
!!                      integer, intent (in) :: blockList (:), 
!!                      real,    intent (in) :: timeStep,
!!                      real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates rays and places them in their initial blocks. This routine calls the
!!  appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : current timestep value
!!  timeSimulation : current simulation time
!!
!!***

subroutine ed_createRays (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ed_createRays
