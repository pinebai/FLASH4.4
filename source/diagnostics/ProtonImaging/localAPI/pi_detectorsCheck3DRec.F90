!!****if* source/diagnostics/ProtonImaging/localAPI/pi_detectorsCheck3DRec
!!
!! NAME
!!
!!  pi_detectorsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call pi_detectorsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected detectors data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!         1) Check, if all detector screen areas are completely outside of the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_detectorsCheck3DRec ()

  implicit none

  return
end subroutine pi_detectorsCheck3DRec
