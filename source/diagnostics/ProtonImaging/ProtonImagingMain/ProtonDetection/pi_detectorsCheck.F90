!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_detectorsCheck
!!
!! NAME
!!
!!  pi_detectorsCheck
!!
!! SYNOPSIS
!!
!!  call pi_detectorsCheck ()
!!
!! DESCRIPTION
!!
!!  Checks the collected detectors data for specific grid geometries. This routine calls the
!!  appropriate subroutines according to the domain geometry specified. If the geometry does not
!!  figure inside this routine it simply does nothing and exits.
!!
!! ARGUMENTS
!!
!!***

subroutine pi_detectorsCheck ()

  use pi_interface,        ONLY : pi_detectorsCheck3DRec

  use ProtonImaging_data,  ONLY : pi_gridGeometry, &
                                  pi_3Din2D

  implicit none

#include "ProtonImaging.h"
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (pi_gridGeometry)

    case (GRID_2DCYLINDRICAL)

      if (pi_3Din2D) then
!          call pi_detectorsCheck2DCyl3D ()
      end if

    case (GRID_3DCARTESIAN)

      call pi_detectorsCheck3DRec ()

    case default

      return

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine pi_detectorsCheck
