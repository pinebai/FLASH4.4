!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleCenterOfExpansion
!!
!! NAME
!!
!!  gr_mpoleCenterOfExpansion
!!
!! SYNOPSIS
!!
!!  gr_mpoleCenterOfExpansion (integer, intent(in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Computes all data related to the center of expansion for the multipoles.
!!  This routine is just the wrapper to call the appropriate routine according
!!  to the geometry present.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

subroutine gr_mpoleCenterOfExpansion (idensvar)

  use gr_mpoleData,      ONLY : gr_mpoleGeometry

  use gr_mpoleInterface, ONLY : gr_mpoleCen3Dcartesian,   &
                                gr_mpoleCen3Dcylindrical, &
                                gr_mpoleCen2Dcylindrical, &
                                gr_mpoleCen2Dspherical,   &
                                gr_mpoleCen1Dspherical
  implicit none

#include "gr_mpole.h"
  
  integer, intent (in) :: idensvar
!  
!
!     ...Select the appropriate subroutine.
!
!
  select case (gr_mpoleGeometry)

    case (GRID_3DCARTESIAN)

          call gr_mpoleCen3Dcartesian   (idensvar)

    case (GRID_3DCYLINDRICAL)

          call gr_mpoleCen3Dcylindrical (idensvar)

    case (GRID_2DCYLINDRICAL)

          call gr_mpoleCen2Dcylindrical (idensvar)

    case (GRID_2DSPHERICAL)

          call gr_mpoleCen2Dspherical   (idensvar)

    case (GRID_1DSPHERICAL)

          call gr_mpoleCen1Dspherical   (idensvar)

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine gr_mpoleCenterOfExpansion
