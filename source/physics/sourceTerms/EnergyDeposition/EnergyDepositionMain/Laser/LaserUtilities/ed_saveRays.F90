!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserUtilities/ed_saveRays
!!
!! NAME
!!
!!  ed_saveRays
!!
!! SYNOPSIS
!!
!!  call ed_saveRays ()
!!
!! DESCRIPTION
!!
!!  This routine saves information of those rays which exited the domain, if requested by
!!  the user. If the user specified no saving action, the routine simply sets all those
!!  rays that exited the domain to nonexistent rays and returns. Otherwise the info is
!!  copied into the saved rays array.
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!***

subroutine ed_saveRays ()

  use Driver_interface,         ONLY : Driver_abortFlash

  use EnergyDeposition_data,    ONLY : ed_maxRayCount,         &
                                       ed_numberOfSavedRays,   &
                                       ed_rayCount,            &
                                       ed_rays,                &
                                       ed_raysSaved,           &
                                       ed_saveOutOfDomainRays

  implicit none

#include "constants.h"
#include "EnergyDeposition.h"

  integer :: blockID
  integer :: ray, rayTag

  real    :: rayPower
  real    :: rayX, rayY, rayZ
!
!
!     ...Save the ray data.
!
!
  if (ed_saveOutOfDomainRays) then

      do ray = 1,ed_rayCount

         blockID = ed_rays (RAY_BLCK,ray)

         if (blockID == RAY_OUTDOMAIN) then

             rayTag   = int (ed_rays (RAY_TAGS,ray))
             rayX     =      ed_rays (RAY_POSX,ray)
             rayY     =      ed_rays (RAY_POSY,ray)
             rayZ     =      ed_rays (RAY_POSZ,ray)
             rayPower =      ed_rays (RAY_POWR,ray)

             ed_numberOfSavedRays = ed_numberOfSavedRays + 1

             if (ed_numberOfSavedRays > ed_maxRayCount) then
                 call Driver_abortFlash ("ed_saveRays: Not enough storage for saved rays array")
             end if

             ed_raysSaved (ed_numberOfSavedRays) % rayX     = rayX
             ed_raysSaved (ed_numberOfSavedRays) % rayY     = rayY
             ed_raysSaved (ed_numberOfSavedRays) % rayZ     = rayZ
             ed_raysSaved (ed_numberOfSavedRays) % rayPower = rayPower
             ed_raysSaved (ed_numberOfSavedRays) % rayTag   = rayTag

             ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)

         end if

      end do
!
!
!     ...Don't save, but mark the rays as nonexistent.
!
!
  else

      do ray = 1,ed_rayCount

         blockID = ed_rays (RAY_BLCK,ray)

         if (blockID == RAY_OUTDOMAIN) then
             ed_rays (RAY_BLCK,ray) = real (NONEXISTENT)
         end if

      end do

  end if
!
!
!     ...Ready!
!
!
  return
end subroutine ed_saveRays
