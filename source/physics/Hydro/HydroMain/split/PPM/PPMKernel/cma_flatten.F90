!!****if* source/physics/Hydro/HydroMain/split/PPM/PPMKernel/cma_flatten
!!
!! NAME
!!
!!  cma_flatten
!!
!! SYNOPSIS
!!
!!  cma_flatten(integer(IN) :: numIntCells,
!!              integer(IN) :: numCells,
!!              integer(IN) :: guard,
!!              real(IN)    :: xn(numCells, hy_numXn),
!!              real(INOUT) :: xnl(numCells, hy_numXn),
!!              real(INOUT) :: xnr(numCells, hy_numXn),
!!              real(INOUT) :: dxn(numCells, hy_numXn),
!!              real(INOUT) :: xn6(numCells, hy_numXn))
!!
!! DESCRIPTION
!!
!!  Look for local extrema of the mass fractions, and apply flattening
!!  to the interpolation profiles in the zones bordering the extrema.
!!  The go over ALL of the abundances, a fix up the interface values
!!  in a manner that should restore the property that they sum to 1.
!!  
!!  Upon exit, the interface values, xnl and xnr may be modified.  
!!  Furthermore, dxn, xn6 will have been updated to be consistent.
!!
!!  For details see "the CMA paper": Plewa & Mueller, 1999 A&A, 342, 179.
!!
!! ARGUMENTS
!!
!! numIntCells :
!! numCells :
!! guard :
!! xn :
!! xnl :
!! xnr :
!! dxn :
!! xn6 :
!!
!!***


subroutine cma_flatten(numIntCells, numCells,guard, xn, xnl, xnr, dxn, xn6)
  use Hydro_data, ONLY : hy_numXn
  implicit none

#include "Flash.h"

  integer, intent(IN) :: numIntCells, numCells, guard
  real, intent(IN),    DIMENSION(numCells, hy_numXn) :: xn
  real, intent(INOUT), DIMENSION(numCells, hy_numXn) :: xnl, xnr, dxn, xn6

  integer :: i, n

  real :: tmp, w_ln, w_rn

  real, DIMENSION(numCells) :: extrema, flatten
  real, DIMENSION(numCells) :: s_L_plus, s_L_minus, s_R_plus, s_R_minus
  real, DIMENSION(numCells) :: delta_L_min, delta_L_max, delta_R_min, delta_R_max
  real, DIMENSION(numCells) :: w_L, w_R

  real, DIMENSION(numCells, hy_numXn) :: sgn_L, sgn_R

! start by looking for extrema, one abundance at a time

  do n = 1, NSPECIES

     do i = guard-1, guard+numIntCells+2
        if ( (xn(i+1,n) - xn(i,n))*(xn(i,n) - xn(i-1,n)) < 0.e0) then
           extrema(i) = 1.0
        else
           extrema(i) = 0.0
        endif
     enddo

! If zone i has an extrema, we want to flatten its neighbors.  Furthermore,
! we only want to flatten the interpolation profiles for zone i if there
! is an extrema in one of its neighbors and NOT in zone i too.

     do i = guard, guard+numIntCells+1
        flatten(i) = 0.5*max(extrema(i-1), 2.0*extrema(i), extrema(i+1))
        flatten(i) = min(1.0, flatten(i))
     enddo

! Now flatten(i) is > 0 if we are to apply flattening in zone i.  It is
! 0.5 if we are to maximally flatten (since one of it's neighboring zones
! contains an extrema), and it is 1.0 if we are to not flatten, since zone
! i itself contains an extrema.

! Apply Eq. 14

     do i = guard, guard+numIntCells+1
        if (flatten(i) > 0.0) then
           tmp = 1.0 - flatten(i)
           
           xnl(i,n) = flatten(i)*xn(i,n) + tmp*xnl(i,n)
           xnr(i,n) = flatten(i)*xn(i,n) + tmp*xnr(i,n)
           
           dxn(i,n) = xnr(i,n) - xnl(i,n)
           xn6(i,n) = 6.0*xn(i,n) - 3.0*(xnl(i,n) + xnr(i,n))
        endif
     enddo

  enddo

! Now we perform the second modification described in the CMA paper.
! Here, we use a variable flattening coefficient to flatten the interface
! values of all of the abundances


! Start by computing two sums, at each interface, Eq. 15.
  do i = guard, guard+numIntCells+1
     s_L_plus(i) = 0.0
     s_L_minus(i) = 0.0

     s_R_plus(i) = 0.0
     s_R_minus(i) = 0.0
  enddo

  do n = 1, NSPECIES
     do i = guard, guard+numIntCells+1
        s_L_plus(i)  = s_L_plus(i)  + max(0.0, xnl(i,n) - xn(i,n))
        s_L_minus(i) = s_L_minus(i) + max(0.0, xn(i,n) - xnl(i,n))
        
        s_R_plus(i)  = s_R_plus(i)  + max(0.0, xnr(i,n) - xn(i,n))
        s_R_minus(i) = s_R_minus(i) + max(0.0, xn(i,n) - xnr(i,n))
     enddo
  enddo

  do i = guard, guard+numIntCells+1
     delta_L_min(i) = min(s_L_plus(i), s_L_minus(i))
     delta_L_max(i) = max(s_L_plus(i), s_L_minus(i))
  
     delta_R_min(i) = min(s_R_plus(i), s_R_minus(i))
     delta_R_max(i) = max(s_R_plus(i), s_R_minus(i))
  enddo

! compute the sign term that is in the w# calculation.  This is different
! for each abundance
  do n = 1, NSPECIES
     do i = guard, guard+numIntCells+1
        sgn_L(i,n) = 0.5*abs(sign(1.0, xnr(i,n) - xnl(i,n)) - &
                             sign(1.0,s_L_plus(i) - s_L_minus(i)))

        sgn_R(i,n) = 0.5*abs(sign(1.0, xnr(i,n) - xnl(i,n)) + &
                             sign(1.0,s_R_plus(i) - s_R_minus(i)))

     enddo
  enddo

! now, fix up the interface values
  do i = guard, guard+numIntCells+1
     w_L(i) = max(0.0, min(1.0, 0.25*(delta_L_max(i) - delta_L_min(i))/ &
                                      delta_L_min(i)))
     w_R(i) = max(0.0, min(1.0, 0.25*(delta_R_max(i) - delta_R_min(i))/ &
                                      delta_R_min(i)))
  enddo

  do n = 1, NSPECIES
     do i = guard, guard+numIntCells+1
        w_Ln = sgn_L(i,n)*w_L(i)
        w_Rn = sgn_R(i,n)*w_R(i)

        xnl(i,n) = w_Ln*xn(i,n) + (1.0 - w_Ln)*xnl(i,n)
        xnr(i,n) = w_Rn*xn(i,n) + (1.0 - w_Rn)*xnr(i,n)
        
        dxn(i,n) = xnr(i,n) - xnl(i,n)
        xn6(i,n) = 6.0*xn(i,n) - 3.0*(xnl(i,n) + xnr(i,n))
     enddo
  enddo

  return
end subroutine cma_flatten

