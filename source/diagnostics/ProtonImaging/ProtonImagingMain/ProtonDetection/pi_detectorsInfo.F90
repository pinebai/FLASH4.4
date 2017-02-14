!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonDetection/pi_detectorsInfo
!!
!! NAME
!!
!!  pi_detectorsInfo
!!
!! SYNOPSIS
!!
!!  call pi_detectorsInfo ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton detectors for specific grid geometries. In here
!!  all detector information is generated that can be generated at initialization. This routine
!!  calls the appropriate subroutines according to the domain geometry specified. If the
!!  geometry does not figure inside this routine it simply does nothing and exits.
!!
!! ARGUMENTS
!!
!!***

subroutine pi_detectorsInfo ()

  use pi_interface,        ONLY : pi_detectorsInfo3DRec

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
          call pi_detectorsInfo3DRec ()
      end if

    case (GRID_3DCARTESIAN)

      call pi_detectorsInfo3DRec ()

    case default

      return

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine pi_detectorsInfo
