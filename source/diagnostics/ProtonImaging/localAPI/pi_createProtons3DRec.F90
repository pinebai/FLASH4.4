!!****if* source/diagnostics/ProtonImaging/localAPI/pi_createProtons3DRec
!!
!! NAME
!!
!!  pi_createProtons3DRec
!!
!! SYNOPSIS
!!
!!  call pi_createProtons3DRec (integer, intent (in) :: blockCount,
!!                              integer, intent (in) :: blockList (:),
!!                              real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates protons and places them in their initial blocks for those geometries consisting
!!  formally of 3D rectangular grids (cartesian). On exit, all protons hitting the domain boundary
!!  have been generated for the current processor. Their block ID's are not ordered as the
!!  outer loop is over all beams.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeSimulation : current simulation time
!!
!! NOTES
!!
!!***

subroutine pi_createProtons3DRec (blockCount, blockList, timeSimulation)

  implicit none

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation

  return
end subroutine pi_createProtons3DRec
