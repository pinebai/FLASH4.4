!!****if* source/diagnostics/ProtonImaging/ProtonImagingMain/ProtonBeams/pi_beamsInfo
!!
!! NAME
!!
!!  pi_beamsInfo
!!
!! SYNOPSIS
!!
!!  call pi_beamsInfo ()
!!
!! DESCRIPTION
!!
!!  Generates information about the proton beams for specific grid geometries. In here
!!  all beam information is generated that can be generated at initialization. This routine
!!  calls the appropriate subroutines according to the domain geometry specified. If the
!!  geometry does not figure inside this routine it simply does nothing and exits.
!!
!! ARGUMENTS
!!
!!***

subroutine pi_beamsInfo ()

  use pi_interface,        ONLY : pi_beamsInfo3DRec

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
          call pi_beamsInfo3DRec ()
      end if

    case (GRID_3DCARTESIAN)

      call pi_beamsInfo3DRec ()

    case default

      return

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine pi_beamsInfo
