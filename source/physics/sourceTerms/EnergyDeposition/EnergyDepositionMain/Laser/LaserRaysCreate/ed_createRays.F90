!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserRaysCreate/ed_createRays
!!
!! NAME
!!
!!  ed_createRays
!!
!! SYNOPSIS
!!
!!  call ed_createRays (integer, intent (in) :: blockCount, 
!!                      integer, intent (in) :: blockList (:), 
!!                      real,    intent (in) :: timeStep,
!!                      real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates rays and places them in their initial blocks. This routine calls the
!!  appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeStep       : current timestep value
!!  timeSimulation : current simulation time
!!
!!***

subroutine ed_createRays (blockCount, blockList, timeStep, timeSimulation)

  use Driver_interface,       ONLY : Driver_abortFlash

  use ed_interface,           ONLY : ed_createRays1DRec,   &
                                     ed_createRays2DCyl3D, &
                                     ed_createRays2DRec,   &
                                     ed_createRays3DRec,   &
                                     ed_initializeRays

  use EnergyDeposition_data,  ONLY : ed_gridGeometry,        &
                                     ed_laser3Din2D,         &
                                     ed_rays,                &
                                     ed_raysMovedIntoDomain

  implicit none

#include "EnergyDeposition.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeStep
  real,    intent (in) :: timeSimulation
!
!
!     ...Initialize the rays.
!
!
  call ed_initializeRays ()
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (ed_gridGeometry)

    case (GRID_1DCARTESIAN)

      call ed_createRays1DRec (blockCount, blockList, timeStep, timeSimulation)

    case (GRID_2DCARTESIAN)

      call ed_createRays2DRec (blockCount, blockList, timeStep, timeSimulation)

    case (GRID_2DCYLINDRICAL)

      if (ed_laser3Din2D) then
          call ed_createRays2DCyl3D (blockCount, blockList, timeStep, timeSimulation)
      else
          call ed_createRays2DRec (blockCount, blockList, timeStep, timeSimulation)
      end if

    case (GRID_3DCARTESIAN)

      call ed_createRays3DRec (blockCount, blockList, timeStep, timeSimulation)

    case default

      call Driver_abortFlash ('[ed_createRays] ERROR: unknown geometry')

  end select
!
!
!     ...Inform the program, that currently all rays are sitting on the domain
!        boundary and have not moved into the domain yet.
!
!
  ed_raysMovedIntoDomain = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine ed_createRays
