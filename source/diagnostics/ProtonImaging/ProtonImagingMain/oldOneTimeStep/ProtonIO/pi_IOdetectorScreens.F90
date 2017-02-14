!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonIO/pi_IOdetectorScreens
!!
!! NAME
!!
!!  pi_IOdetectorScreens
!!
!! SYNOPSIS
!!
!!  call pi_IOdetectorScreens ()
!!
!! DESCRIPTION
!!
!!  This routine stores IO protons that travel along the edges of each detector
!!  screen into the corresponding IO proton arrays. Five x,y,z positions need to
!!  be recorded in order for each detector IO proton to complete the screen
!!  perimeter. Each such detector IO proton is going to be assigned a tag that
!!  exceeds the current maximum tag value.
!!
!!  Only the master processor generates the detector IO protons.
!!
!!***

subroutine pi_IOdetectorScreens ()
  
  use Driver_interface,    ONLY : Driver_abortFlash

  use ProtonImaging_data,  ONLY : pi_detectors,                &
                                  pi_globalMe,                 &
                                  pi_IOmaxPointsPerBlock,      &
                                  pi_IOmaxProtonCount,         &
                                  pi_IOprotonCount,            &
                                  pi_IOprotonPointCount,       &
                                  pi_IOprotonPoints,           &
                                  pi_IOprotonTags,             &
                                  pi_numberOfDetectors,        &
                                  pi_tagMax

  implicit none

#include "constants.h"
#include "Flash.h"
#include "ProtonImaging.h"

  integer :: detector
!
!
!     ...Return immediately, if the detector IO proton points cannot be stored.
!
!
  if (pi_IOmaxPointsPerBlock < 5) then
      print *, '[pi_IOdetectorScreens] The detector screen IO protons need 5 points, &
                 but have only ', pi_IOmaxPointsPerBlock, ' in buffer space'
      return
  end if

  if (pi_IOmaxProtonCount < pi_numberOfDetectors) then
      print *, '[pi_IOdetectorScreens] # of detectors is ', pi_numberOfDetectors, &
               ' but the max dimension for IO proton count is ',pi_IOmaxProtonCount, &
               ' in buffer space'
      return
  end if
!
!
!     ...Loop over all detectors. Perform the detector edge 'travel' by starting
!        at the upper right corner, followed by the lower right, lower left, upper
!        left corners and ending again at the upper right corner.
!
!
  pi_IOprotonCount = 0

  if (pi_globalMe /= MASTER_PE) then
      return
  end if

  do detector = 1, pi_numberOfDetectors

     pi_tagMax = pi_tagMax + 1
     pi_IOprotonCount = pi_IOprotonCount + 1

     pi_IOprotonTags       (   pi_IOprotonCount       ) = pi_tagMax
     pi_IOprotonPointCount (   pi_IOprotonCount       ) = 5
     pi_IOprotonPoints     (1, pi_IOprotonCount, IAXIS) = pi_detectors (detector) % cornerUpperRightX
     pi_IOprotonPoints     (1, pi_IOprotonCount, JAXIS) = pi_detectors (detector) % cornerUpperRightY
     pi_IOprotonPoints     (1, pi_IOprotonCount, KAXIS) = pi_detectors (detector) % cornerUpperRightZ
     pi_IOprotonPoints     (2, pi_IOprotonCount, IAXIS) = pi_detectors (detector) % cornerLowerRightX
     pi_IOprotonPoints     (2, pi_IOprotonCount, JAXIS) = pi_detectors (detector) % cornerLowerRightY
     pi_IOprotonPoints     (2, pi_IOprotonCount, KAXIS) = pi_detectors (detector) % cornerLowerRightZ
     pi_IOprotonPoints     (3, pi_IOprotonCount, IAXIS) = pi_detectors (detector) % cornerLowerLeftX
     pi_IOprotonPoints     (3, pi_IOprotonCount, JAXIS) = pi_detectors (detector) % cornerLowerLeftY
     pi_IOprotonPoints     (3, pi_IOprotonCount, KAXIS) = pi_detectors (detector) % cornerLowerLeftZ
     pi_IOprotonPoints     (4, pi_IOprotonCount, IAXIS) = pi_detectors (detector) % cornerUpperLeftX
     pi_IOprotonPoints     (4, pi_IOprotonCount, JAXIS) = pi_detectors (detector) % cornerUpperLeftY
     pi_IOprotonPoints     (4, pi_IOprotonCount, KAXIS) = pi_detectors (detector) % cornerUpperLeftZ
     pi_IOprotonPoints     (5, pi_IOprotonCount, IAXIS) = pi_detectors (detector) % cornerUpperRightX
     pi_IOprotonPoints     (5, pi_IOprotonCount, JAXIS) = pi_detectors (detector) % cornerUpperRightY
     pi_IOprotonPoints     (5, pi_IOprotonCount, KAXIS) = pi_detectors (detector) % cornerUpperRightZ
     
  end do
!
!
!    ...Ready!
!
!
  return
end subroutine pi_IOdetectorScreens
