!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserIO/ed_laserIOInit
!!
!! NAME
!!
!!  ed_laserIOInit
!!
!! SYNOPSIS
!!
!!  call ed_laserIOInit ()
!!
!! DESCRIPTION
!!
!!  This routine initializes variables for outputing laser rays.
!!
!!***

subroutine ed_laserIOInit ()
  
  use Driver_interface,             ONLY : Driver_abortFlash

  use EnergyDeposition_data,        ONLY : ed_beams,                       &
                                           ed_laserIOMaxNumberOfPositions, &
                                           ed_laserIOMaxNumberOfRays,      &
                                           ed_laserIONumberOfPositions,    &
                                           ed_laserIORayFrequency,         &
                                           ed_laserIORayPositions,         &
                                           ed_laserIORayPower,             &
                                           ed_laserIORayTags,              &
                                           ed_numberOfBeams,               &
                                           ed_useLaserIO

  use RuntimeParameters_interface,  ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"

  integer :: maxNumberIORays
  integer :: totalRays
  
  call RuntimeParameters_get ('ed_useLaserIO',     ed_useLaserIO)

  if(.not. ed_useLaserIO) return
!
!
!     ...Get the needed external data.
!
!
  call RuntimeParameters_get ("ed_laserIOMaxNumberOfPositions",  ed_laserIOMaxNumberOfPositions)
  call RuntimeParameters_get ("ed_laserIOMaxNumberOfRays",       ed_laserIOMaxNumberOfRays     )

  if (ed_laserIOMaxNumberOfPositions < 0) then
      call Driver_abortFlash ("Must set ed_laserIOMaxNumberOfPositions")
  end if

  if (ed_laserIOMaxNumberOfRays < 0) then
      call Driver_abortFlash ("Must set ed_laserIOMaxNumberOfRays")
  end if
!
!
!     ...Sum over all beams to determine the total number of rays that will be launched.
!        This will determine the ray write out frequency. Adjust the maximum number of
!        IO rays, such that they correspond with the ray write out frequency.
!
!
  totalRays = sum (ed_beams (1:ed_numberOfBeams) % numberOfRays)

  ed_laserIORayFrequency = totalRays / ed_laserIOMaxNumberOfRays
  ed_laserIORayFrequency = max (1,ed_laserIORayFrequency)

  maxNumberIORays = totalRays / ed_laserIORayFrequency
!
!
!     ...Allocate the needed IO arrays.
!
!        The actual number of rays that will be stored is equal to:
!
!              numberIORays = totalRays / ed_laserIORayFrequency
!
!        Note that this WILL break if additional rays are launched on a
!        time step (i.e. if more than totalRays rays are launched). In this
!        case, the arrays that are allocated below will not be big
!        enough. totalRays must represent the maximum number of rays that
!        will be launched on a step.
!
!
  allocate (ed_laserIORayTags           (1:maxNumberIORays                                        ))
  allocate (ed_laserIONumberOfPositions (1:maxNumberIORays                                        ))
  allocate (ed_laserIORayPower          (1:maxNumberIORays, 1:ed_laserIOMaxNumberOfPositions      ))
  allocate (ed_laserIORayPositions      (1:maxNumberIORays, 1:ed_laserIOMaxNumberOfPositions, MDIM))
!
!
!    ...Ready!
!
!
  return
end subroutine ed_laserIOInit
