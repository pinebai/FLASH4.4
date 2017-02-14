!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_detectorsCheck3DRec
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

  use Driver_interface,      ONLY : Driver_abortFlash

  use ProtonImaging_data,    ONLY : pi_detectors,         &
                                    pi_numberOfDetectors, &
                                    pi_xminDomain,        &
                                    pi_xmaxDomain,        &
                                    pi_yminDomain,        &
                                    pi_ymaxDomain,        &
                                    pi_zminDomain,        &
                                    pi_zmaxDomain
  
  implicit none

#include "Flash.h"
#include "ProtonImaging.h"
#include "constants.h"

  logical :: inDomain

  integer :: detector

  real    :: cornerXmin, cornerXmax
  real    :: cornerYmin, cornerYmax
  real    :: cornerZmin, cornerZmax
!
!
!     ...Since the detector screens are square, the extremum values in each direction
!        correspond to one of the four corners. For each X,Y,Z component of the corners
!        we need to find the maximum and the minimum between the four corner values and
!        compare it to the domain extensions.
!
!
  do detector = 1, pi_numberOfDetectors

     cornerXmin = min ( pi_detectors (detector) % cornerUpperRightX, &
                        pi_detectors (detector) % cornerUpperLeftX , &
                        pi_detectors (detector) % cornerLowerRightX, &
                        pi_detectors (detector) % cornerLowerLeftX   )

     cornerXmax = max ( pi_detectors (detector) % cornerUpperRightX, &
                        pi_detectors (detector) % cornerUpperLeftX , &
                        pi_detectors (detector) % cornerLowerRightX, &
                        pi_detectors (detector) % cornerLowerLeftX   )

     cornerYmin = min ( pi_detectors (detector) % cornerUpperRightY, &
                        pi_detectors (detector) % cornerUpperLeftY , &
                        pi_detectors (detector) % cornerLowerRightY, &
                        pi_detectors (detector) % cornerLowerLeftY   )

     cornerYmax = max ( pi_detectors (detector) % cornerUpperRightY, &
                        pi_detectors (detector) % cornerUpperLeftY , &
                        pi_detectors (detector) % cornerLowerRightY, &
                        pi_detectors (detector) % cornerLowerLeftY   )

     cornerZmin = min ( pi_detectors (detector) % cornerUpperRightZ, &
                        pi_detectors (detector) % cornerUpperLeftZ , &
                        pi_detectors (detector) % cornerLowerRightZ, &
                        pi_detectors (detector) % cornerLowerLeftZ   )

     cornerZmax = max ( pi_detectors (detector) % cornerUpperRightZ, &
                        pi_detectors (detector) % cornerUpperLeftZ , &
                        pi_detectors (detector) % cornerLowerRightZ, &
                        pi_detectors (detector) % cornerLowerLeftZ   )
!
!
!     ...Check, if detector screen area touches the domain.
!
!
     inDomain =      (cornerXmin >= pi_xminDomain) &
               .and. (cornerYmin >= pi_yminDomain) &
               .and. (cornerZmin >= pi_zminDomain) &
               .and. (cornerXmax <= pi_xmaxDomain) &
               .and. (cornerYmax <= pi_ymaxDomain) &
               .and. (cornerZmax <= pi_zmaxDomain)

     if (inDomain) then
         call Driver_abortFlash ("pi_detectorsCheck3DRec: Detector screen (partially) inside the domain!")
     end if

  end do
!
!
!     ...Ready!
!
!
  return
end subroutine pi_detectorsCheck3DRec
