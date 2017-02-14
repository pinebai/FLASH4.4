!!****if* source/diagnostics/ProtonImaging/localAPI/pi_beamsCheck3DRec
!!
!! NAME
!!
!!  pi_beamsCheck3DRec
!!
!! SYNOPSIS
!!
!!  call pi_beamsCheck3DRec ()
!!
!! DESCRIPTION
!!
!!  Checks the collected beams data for those geometries consisting formally of 3D rectangular
!!  grids (cartesian). All checks which depend on domain grid details should go in here.
!!  Currently it contains the following:
!!
!!         1) Check, if all beam circular lens areas are completely outside the domain.
!!
!!  Since the target area of the proton beam is merely used to construct the direction of
!!  the individual protons, there is no need to check if this area is properly located wrt
!!  to the domain.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine pi_beamsCheck3DRec ()

  implicit none

  return
end subroutine pi_beamsCheck3DRec
