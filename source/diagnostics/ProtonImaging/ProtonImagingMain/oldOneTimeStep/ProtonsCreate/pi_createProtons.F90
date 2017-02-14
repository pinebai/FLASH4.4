!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsCreate/pi_createProtons
!!
!! NAME
!!
!!  pi_createProtons
!!
!! SYNOPSIS
!!
!!  call pi_createProtons (integer, intent (in) :: blockCount, 
!!                         integer, intent (in) :: blockList (:),
!!                         real,    intent (in) :: timeSimulation)
!!
!! DESCRIPTION
!!
!!  Generates protons and places them in their initial blocks. This routine calls the
!!  appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!  blockCount     : Number of blocks on current processor
!!  blockList      : All block ID numbers
!!  timeSimulation : current simulation time
!!
!!***

subroutine pi_createProtons (blockCount, blockList, timeSimulation)

  use Driver_interface,    ONLY : Driver_abortFlash

  use pi_interface,        ONLY : pi_createProtons3DRec

  use ProtonImaging_data,  ONLY : pi_gridGeometry,           &
                                  pi_3Din2D,                 &
                                  pi_protonCount,            &
                                  pi_protons,                &
                                  pi_protonsMovedIntoDomain

  implicit none

#include "ProtonImaging.h"

  integer, intent (in) :: blockCount
  integer, intent (in) :: blockList (1:blockCount)
  real,    intent (in) :: timeSimulation
!
!
!     ...Initialize the protons count.
!
!
  pi_protonCount = 0
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (pi_gridGeometry)

    case (GRID_2DCYLINDRICAL)

      if (pi_3Din2D) then
          call Driver_abortFlash ('[pi_createProtons] ERROR: unsupported/unknown geometry')
!          call pi_createProtons2DCyl3D (blockCount, blockList, timeSimulation)
      end if

    case (GRID_3DCARTESIAN)

      call pi_createProtons3DRec (blockCount, blockList, timeSimulation)

    case default

      call Driver_abortFlash ('[pi_createProtons] ERROR: unsupported/unknown geometry')

  end select
!
!
!     ...Inform the program, that currently all protons are sitting on the domain
!        boundary and have not moved into the domain yet.
!
!
  pi_protonsMovedIntoDomain = .false.
!
!
!    ...Ready!
!
!
  return
end subroutine pi_createProtons
