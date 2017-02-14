!!****if* source/physics/sourceTerms/EnergyDeposition/localAPI/ed_beamsCheck3DRec
!!
!! NAME
!!
!!  ed_beamsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call ed_beamsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!         1) Check, if all beam elliptical target areas are completely within the domain.
!!         2) Check, if all beam elliptical lens areas are completely outside the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_beamsCheck3DRec ()

  implicit none

  return
end subroutine ed_beamsCheck3DRec
