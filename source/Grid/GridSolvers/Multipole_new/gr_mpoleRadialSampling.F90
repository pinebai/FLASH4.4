!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleRadialSampling
!!
!! NAME
!!
!!  gr_mpoleRadialSampling
!!
!! SYNOPSIS
!!
!!  gr_mpoleRadialSampling ()
!!
!! DESCRIPTION
!!
!!  Determines the radial sampling for accumulating the moments.
!!  It defines the radial bins into which each moment falls. If the
!!  radial sampling is chosen unwisely, lots of memory for holding
!!  empty moment bins might be wasted. This routine calls the appropriate
!!  subroutines according to the geometry specified.
!!
!!***

subroutine gr_mpoleRadialSampling ()

  use gr_mpoleInterface, ONLY : gr_mpoleRad3Dcartesian,   &
                                gr_mpoleRad3Dcylindrical, &
                                gr_mpoleRad2Dcylindrical, &
                                gr_mpoleRad2Dspherical,   &
                                gr_mpoleRad1Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleGeometry

  implicit none

#include "gr_mpole.h"
!
!
!     ...Select the appropriate subroutine.
!
!
  select case (gr_mpoleGeometry)

    case (GRID_3DCARTESIAN)

          call gr_mpoleRad3Dcartesian   ()

    case (GRID_3DCYLINDRICAL)

          call gr_mpoleRad3Dcylindrical ()

    case (GRID_2DCYLINDRICAL)

          call gr_mpoleRad2Dcylindrical ()

    case (GRID_2DSPHERICAL)

          call gr_mpoleRad2Dspherical   ()

    case (GRID_1DSPHERICAL)

          call gr_mpoleRad1Dspherical   ()

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpoleRadialSampling
