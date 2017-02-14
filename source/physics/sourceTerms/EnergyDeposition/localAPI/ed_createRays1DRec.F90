!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_createRays1DRec
!!
!! NAME
!!
!!  ed_createRays1DRec
!!
!! SYNOPSIS
!!
!!  call ed_createRays1DRec (integer, intent (in) :: blockCount, 
!!                           integer, intent (in) :: blockList (:), 
!!                           real,    intent (in) :: timeStep,   
!!                           real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates one ray per beam and places it in its initial block for those geometries consisting
!!  formally of 1D rectangular grids (cartesian + spherical). On exit, all rays hitting
!!  the domain boundary have been generated for the current processor. Their block ID's are
!!  not ordered as the outer loop is over all beams.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : current timestep value
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!  Since this routine deals with 1D rectangular type of grids, the following unique features
!!  will occur:
!!
!!   1) The number of rays created must be equal to the number of active beams
!!   2) The number of distinct block ID's can be at most 2, corresponding
!!      to the two outer domain blocks.
!!
!!***

subroutine ed_createRays1DRec (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ed_createRays1DRec
