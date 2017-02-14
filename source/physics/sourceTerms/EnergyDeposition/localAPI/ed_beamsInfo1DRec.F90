!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsInfo1DRec
!!
!! NAME
!!
!!  ed_beamsInfo1DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsInfo1DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the beams for those geometries consisting formally of 1D
!!  rectangular grids (cartesian + spherical). In here all beam information is generated
!!  that can be generated at initialization. Currently it contains the following:
!!
!!         1) Calculate the lens to target distance and incorporate some
!!            additional info (time duration) regarding the pulses associated
!!            with each beam.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsInfo1DRec ()

  implicit none

  return
end subroutine ed_beamsInfo1DRec
