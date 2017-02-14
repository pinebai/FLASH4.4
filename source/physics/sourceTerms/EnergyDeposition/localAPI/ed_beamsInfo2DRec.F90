!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsInfo2DRec
!!
!! NAME
!!
!!  ed_beamsInfo2DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo2DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those geometries consisting formally of 2D
!!  rectangular grids (cartesian + cylindrical). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!!         2) Calculate basic information about the linear beam ray grid.
!!
!!         3) Calculate the unit vectors associated with both the linear axis
!!            in the local target coordinate system (i.e. with origin on
!!            the center of the target line).
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo2DRec ()

  implicit none

  return
end subroutine ed_beamsInfo2DRec
