!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsCheck1DRec
!!
!! NAME
!!
!!  ed_beamsCheck1DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck1DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 1D rectangular
!!  grids (cartesian + spherical). All checks which depend on domain grid details should go
!!  in here. Currently it contains the following:
!!
!!         1) Check, if all beam target points are completely within the domain.
!!         2) Check, if all beam lens points are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck1DRec ()

  implicit none

  return
end subroutine ed_beamsCheck1DRec
