!!****if* source/diagnostics/ProtonImaging/localAPI/pi_beamsInfo3DRec
!!
!! NAME
!!
!!  pi_beamsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call pi_beamsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton beams for those domain geometries consisting
!!  formally of 3D rectangular grids (cartesian). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance.
!!
!!         2) Calculate basic information about the planar circular beam proton grid.
!!
!!         3) Calculate two orthogonal unit vectors within the circular plane of
!!            of the lens located at the lens' center.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!***

subroutine pi_beamsInfo3DRec ()

  implicit none

  return
end subroutine pi_beamsInfo3DRec
