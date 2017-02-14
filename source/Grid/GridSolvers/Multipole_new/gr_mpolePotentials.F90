!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpolePotentials
!!
!! NAME
!!
!!  gr_mpolePotentials
!!
!! SYNOPSIS
!!
!!  gr_mpolePotentials  (integer, intent(in) :: ipotvar,
!!                       real,    intent(in) :: Poisson_factor )
!!
!! DESCRIPTION
!!
!!  Computes the potential field using the mass moments already
!!  calculated. On output tha variable indexed by ipotvar contains
!!  the potential. The calculations are entirely local to each
!!  processor, since each processor has a local copy of the moments.
!!  This routine calls the appropriate subroutines according to
!!  the domain geometry specified.
!!
!! ARGUMENTS
!!
!!  ipotvar        : index to variable containing the potential
!!  Poisson_factor : the factor in front of the Poisson equation
!!
!!***

subroutine gr_mpolePotentials (ipotvar,Poisson_factor)

  use gr_mpoleInterface, ONLY : gr_mpolePot3Dcartesian,   &
                                gr_mpolePot3Dcylindrical, &
                                gr_mpolePot2Dcylindrical, &
                                gr_mpolePot2Dspherical,   &
                                gr_mpolePot1Dspherical

  use gr_mpoleData,      ONLY : gr_mpoleGravityConstant, &
                                gr_mpoleFourPiInv,       &
                                gr_mpoleGeometry

  implicit none

#include "gr_mpole.h"

  integer, intent (in) :: ipotvar
  real,    intent (in) :: Poisson_factor
!  
!
!     ...Calculate the gravitational constant.
!
!
!$omp single
  gr_mpoleGravityConstant = Poisson_factor * gr_mpoleFourPiInv
!$omp end single

  !OpenMP implicit barrier at the end of the single is needed.
!
!
!     ...Select the appropriate subroutine.
!
!

  select case (gr_mpoleGeometry)

     case (GRID_3DCARTESIAN)

           call gr_mpolePot3Dcartesian   (ipotvar)

     case (GRID_3DCYLINDRICAL)

           call gr_mpolePot3Dcylindrical (ipotvar)

     case (GRID_2DCYLINDRICAL)

           call gr_mpolePot2Dcylindrical (ipotvar)

     case (GRID_2DSPHERICAL)

           call gr_mpolePot2Dspherical   (ipotvar)

     case (GRID_1DSPHERICAL)

           call gr_mpolePot1Dspherical   (ipotvar)

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine gr_mpolePotentials
