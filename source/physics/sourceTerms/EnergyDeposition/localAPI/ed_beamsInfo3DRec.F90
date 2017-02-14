!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsInfo3DRec
!!
!! NAME
!!
!!  ed_beamsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those geometries consisting formally of 3D
!!  rectangular grids (cartesian). In here all beam information is generated that can be
!!  generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!!         2) Calculate basic information about the planar beam ray grid.
!!
!!         3) Calculate the unit vectors associated with both the semiaxes
!!            in the local target coordinate system (i.e. with origin on
!!            the center of the ellipse).
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo3DRec ()

  implicit none

  return
end subroutine ed_beamsInfo3DRec
