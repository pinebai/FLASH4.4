!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRayTrace/ed_traceRays
!!
!! NAME
!!
!!  ed_traceRays
!!
!! SYNOPSIS
!!
!!  call ed_traceRays (real, intent (in) :: timeStep)
!!
!! DESCRIPTION
!!
!!  Follows the collection of all rays at this stage through one block and deposits their energy
!!  in those blocks. This routine calls the appropriate subroutines according to the domain
!!  grid geometry specified.
!!
!! ARGUMENTS
!!
!!  timeStep : Current time step
!!
!!***

subroutine ed_traceRays (timeStep) 

  use Driver_interface,       ONLY : Driver_abortFlash

  use ed_interface,           ONLY : ed_traceRays1DRec,   &
                                     ed_traceRays2DCyl3D, &
                                     ed_traceRays2DRec,   &
                                     ed_traceRays3DRec

  use EnergyDeposition_data,  ONLY : ed_gridGeometry,               &
                                     ed_laser3Din2D,                &
                                     ed_laserIONumberOfRaysWritten, &
                                     ed_rayCount,                   &
                                     ed_raysMovedIntoDomain

  implicit none

#include "EnergyDeposition.h"

  real, intent (in) :: timeStep
!
!
!     ...Reset the number of rays to write out to zero.
!
!
  ed_laserIONumberOfRaysWritten = 0
!
!
!     ...Select the appropriate subroutine.
!
!
  if (ed_rayCount > 0) then

      select case (ed_gridGeometry)

      case (GRID_1DCARTESIAN)

        call ed_traceRays1DRec (timeStep)

      case (GRID_2DCARTESIAN)

        call ed_traceRays2DRec (timeStep)

      case (GRID_2DCYLINDRICAL)

        if (ed_laser3Din2D) then
            call ed_traceRays2DCyl3D (timeStep)
        else
            call ed_traceRays2DRec (timeStep)
        end if

      case (GRID_3DCARTESIAN)

        call ed_traceRays3DRec (timeStep)

      case default

        call Driver_abortFlash ('[ed_traceRays] ERROR: unknown geometry')

      end select

  end if
!
!
!     ...After at least one call to the ray tracing routine, the rays have
!        moved into the domain.
!
!
  ed_raysMovedIntoDomain = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine ed_traceRays
