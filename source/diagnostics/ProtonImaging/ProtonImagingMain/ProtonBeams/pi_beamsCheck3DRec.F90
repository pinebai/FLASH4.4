!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonBeams/pi_beamsCheck3DRec
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
!!         1) Check, if all beam capsule spheres are completely outside the domain.
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

  use Driver_interface,      ONLY : Driver_abortFlash

  use ProtonImaging_data,    ONLY : pi_beams,         &
                                    pi_numberOfBeams, &
                                    pi_xminDomain,    &
                                    pi_xmaxDomain,    &
                                    pi_yminDomain,    &
                                    pi_ymaxDomain,    &
                                    pi_zminDomain,    &
                                    pi_zmaxDomain
  
  implicit none

#include "Flash.h"
#include "ProtonImaging.h"
#include "constants.h"

  logical :: outOfDomain

  integer :: beam

  real    :: capsuleRadius
  real    :: capsuleX, capsuleY, capsuleZ
!
!
!     ...Loop over all beams.
!
!
  do beam = 1, pi_numberOfBeams

     capsuleRadius = pi_beams (beam) % capsuleRadius
     capsuleX      = pi_beams (beam) % capsuleX
     capsuleY      = pi_beams (beam) % capsuleY
     capsuleZ      = pi_beams (beam) % capsuleZ
!
!
!     ...Check, if beam capsule sphere is completely outside of the domain.
!
!
     outOfDomain =     (capsuleX - capsuleRadius > pi_xmaxDomain) &
                  .or. (capsuleY - capsuleRadius > pi_ymaxDomain) &
                  .or. (capsuleZ - capsuleRadius > pi_zmaxDomain) &
                  .or. (capsuleX + capsuleRadius < pi_xminDomain) &
                  .or. (capsuleY + capsuleRadius < pi_yminDomain) &
                  .or. (capsuleZ + capsuleRadius < pi_zminDomain)

     if (.not.outOfDomain) then
         call Driver_abortFlash ("pi_beamsCheck3DRec: Beam capsule (partially) inside the domain!")
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_beamsCheck3DRec
