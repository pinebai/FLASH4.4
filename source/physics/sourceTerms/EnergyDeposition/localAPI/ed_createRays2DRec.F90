!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_createRays2DRec
!!
!! NAME
!!
!!  ed_createRays2DRec
!!
!! SYNOPSIS
!!
!!  call ed_createRays2DRec (integer, intent (in) :: blockCount, 
!!                           integer, intent (in) :: blockList (:), 
!!                           real,    intent (in) :: timeStep,   
!!                           real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates rays and places them in their initial blocks for those geometries consisting
!!  formally of 2D rectangular grids (cartesian + cylindrical). On exit, all rays hitting
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
!!***

subroutine ed_createRays2DRec (blockCount, blockList, timeStep, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation

  return
end subroutine ed_createRays2DRec
