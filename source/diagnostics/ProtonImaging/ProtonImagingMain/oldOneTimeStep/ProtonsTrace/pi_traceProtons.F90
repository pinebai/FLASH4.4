!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/oldOneTimeStep/ProtonsTrace/pi_traceProtons
!!
!! NAME
!!
!!  pi_traceProtons
!!
!! SYNOPSIS
!!
!!  call pi_traceProtons ()
!!
!! DESCRIPTION
!!
!!  Follows the collection of all protons at this stage through one block. This routine calls
!!  the appropriate subroutines according to the domain grid geometry specified.
!!
!! ARGUMENTS
!!
!!***

subroutine pi_traceProtons () 

  use Driver_interface,    ONLY : Driver_abortFlash

  use pi_interface,        ONLY : pi_traceProtons3DRec

  use ProtonImaging_data,  ONLY : pi_gridGeometry,             &
                                  pi_3Din2D,                   &
                                  pi_IOprotonCount,            &
                                  pi_protonCount,              &
                                  pi_protonsMovedIntoDomain

  implicit none

#include "ProtonImaging.h"
!
!
!     ...Reset the number of IO protons stored per block to zero.
!
!
  pi_IOprotonCount = 0
!
!
!     ...Select the appropriate subroutine.
!
!
  if (pi_protonCount > 0) then

      select case (pi_gridGeometry)

      case (GRID_2DCYLINDRICAL)

        if (pi_3Din2D) then
            call Driver_abortFlash ('[pi_traceProtons] ERROR: no pi_traceProtons2DCyl3D routine')
        end if

      case (GRID_3DCARTESIAN)

        call pi_traceProtons3DRec ()

      case default

        call Driver_abortFlash ('[pi_traceProtons] ERROR: unsupported/unknown geometry')

      end select

  end if
!
!
!     ...After at least one call to the proton tracing routine, the protons have
!        moved into the domain.
!
!
  pi_protonsMovedIntoDomain = .true.
!
!
!    ...Ready!
!
!
  return
end subroutine pi_traceProtons
