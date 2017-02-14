!!****f* source/physics/sourceTerms/EnergyDeposition/EnergyDeposition
!!
!! NAME
!!
!!  EnergyDeposition
!!
!! SYNOPSIS
!!
!!  call EnergyDeposition (integer, intent (in)           :: blockCount, 
!!                         integer, intent (in)           :: blockList (:), 
!!                         real,    intent (in)           :: timeStep,
!!                         real,    intent (in)           :: timeSimulation,
!!                         integer, intent (in), optional :: passSplitDriver)
!!
!! DESCRIPTION
!!
!!  Compute the energy deposited due to irradiation by a laser.
!!  This is the driver routine for computing the energy deposition
!!  during one timestep. It is assumed that the domain is of
!!  dimensions no larger than it takes for light to travel during
!!  this timestep. Hence we follow all rays through the complete
!!  domain during the timestep.
!!
!!  The code operates in a loop where a local routine computes the 
!!  energy deposition and traverses the rays through individual blocks. 
!!  The loop terminates when all the rays have been processed.
!!
!! ARGUMENTS
!!
!!  blockCount      : Number of blocks on current processor
!!  blockList       : All block ID numbers
!!  timeStep        : current timestep value
!!  timeSimulation  : current simulation time
!!  passSplitDriver : indicates first/second half of time step for split driver
!!
!! NOTES
!!          
!!  The current implementation assumes presence of one or more laser 
!!  beams and their path is computed using the geometric optics assumptions.
!!
!!***

subroutine EnergyDeposition (blockCount, blockList, timeStep, timeSimulation, passSplitDriver)

  implicit none

  integer, intent (in)           :: blockCount
  integer, intent (in)           :: blockList (1:blockCount)
  real,    intent (in)           :: timeStep
  real,    intent (in)           :: timeSimulation
  integer, intent (in), optional :: passSplitDriver

  return
end subroutine EnergyDeposition
