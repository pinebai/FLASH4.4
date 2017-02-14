!!****if* source/diagnostics/ProtonImaging/localAPI/pi_detectorsInfo3DRec
!!
!! NAME
!!
!!  pi_detectorsInfo3DRec
!!
!! SYNOPSIS
!!
!!  call pi_detectorsInfo3DRec ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton detectors placed in 3D rectangular (cartesian)
!!  space. In here all detector information is generated that can be generated at initialization.
!!  Currently it contains the following:
!!
!!         1) Calculate 1/2 and inverse of side length (for efficiency purposes).
!!
!!         2) Calculate the normal vector (if detector is aligned wrt to a proton beam).
!!
!!         3) Calculate the unit normal vector.
!!
!!         4) Calculate two orthogonal unit vectors within the square screen plane of
!!            each detector, which will serve as a 2D cartesian local basis for the
!!            screen, located at the screen center.
!!
!! ARGUMENTS
!!
!!  none
!!
!! NOTES
!!
!!***

subroutine pi_detectorsInfo3DRec ()

  implicit none

  return
end subroutine pi_detectorsInfo3DRec
