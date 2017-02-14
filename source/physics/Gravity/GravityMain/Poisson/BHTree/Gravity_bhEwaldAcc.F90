!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhEwaldAcc
!!
!! NAME
!!
!!  Gravity_bhEwaldAcc
!!
!! 
!! SYNOPSIS
!!
!!   real field = Gravity_bhEwaldAcc(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z,
!!                           real(in)    :: drAbsInv
!!        )
!!
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field for acceleration and returns its 
!!   value for the point x,y,z corrected to x_i/r**3 term
!!
!! ARGUMENTS
!!
!!  x   : x-coordinate of the point where the Ewald field is determined
!!  y   : y-coordinate of the point where the Ewald field is determined
!!  z   : z-coordinate of the point where the Ewald field is determined
!!  drAbsInv   :  1.0/sqrt(x**2+y**2+z**2)
!!
!! RESULT
!!
!!  Value of the normalised acceleration in a given point.
!!
!! NOTES
!!
!!***

function Gravity_bhEwaldAcc(x, y, z, drAbsInv)

  use Gravity_data, ONLY : grv_bhDxI, &
      grv_bhTreeEwald, grv_L2inv, grv_margin, &
      grv_bhEwald_periodicity, grv_bhLayerAccConst !, nin, nout
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  real, intent(IN) :: x, y, z, drAbsInv

  integer :: i, j, k
  real :: b, drAbsInv2, drAbsInv3, dradInv2
  real :: Gravity_bhEwaldAcc(IAXIS:KAXIS)


  i = floor(0.5D0+x*grv_bhDxI)
  j = floor(0.5D0+y*grv_bhDxI)
  k = floor(0.5D0+z*grv_bhDxI)

  drAbsInv2 = drAbsInv**2
! if the point is outside the Ewald field, analytical approximation is applied
  if ((max(i,j,k)).GT.grv_margin) then
    select case (grv_bhEwald_periodicity)
      case (1)
! Taylor expansion of 1/R**2 for x/R << 1, R is now radial distance to the cylinder axis
        b = (drAbsInv*x)**2
        dradInv2 = grv_L2inv*drAbsInv2*(1.0+b*(1.0+b*(1.0+b)))
        Gravity_bhEwaldAcc(1) = 0.0
        Gravity_bhEwaldAcc(2) = y*dradInv2
        Gravity_bhEwaldAcc(3) = z*dradInv2
      case (2)
        b = (drAbsInv*y)**2
        dradInv2 = grv_L2inv*drAbsInv2*(1.0+b*(1.0+b*(1.0+b)))
        Gravity_bhEwaldAcc(1) = x*dradInv2
        Gravity_bhEwaldAcc(2) = 0.0D0
        Gravity_bhEwaldAcc(3) = z*dradInv2
      case (4)
        b = (drAbsInv*z)**2
        dradInv2 = grv_L2inv*drAbsInv2*(1.0+b*(1.0+b*(1.0+b)))
        Gravity_bhEwaldAcc(1) = x*dradInv2
        Gravity_bhEwaldAcc(2) = y*dradInv2
        Gravity_bhEwaldAcc(3) = 0.0D0
      case (3)
        Gravity_bhEwaldAcc(1) = 0.0
        Gravity_bhEwaldAcc(2) = 0.0
!        Gravity_bhEwaldAcc(3) = -grv_bhLayerAccConst
        Gravity_bhEwaldAcc(3) = grv_bhLayerAccConst
      case (5)
        Gravity_bhEwaldAcc(1) = 0.0
        Gravity_bhEwaldAcc(2) = grv_bhLayerAccConst
        Gravity_bhEwaldAcc(3) = 0.0
      case (6)
        Gravity_bhEwaldAcc(1) = grv_bhLayerAccConst
        Gravity_bhEwaldAcc(2) = 0.0
        Gravity_bhEwaldAcc(3) = 0.0
      case (7)
        call Driver_abortFlash("This place must not be accessed")
    end select

  else
! interpolation and adding the x/r**3 term

    drAbsInv3 = drAbsInv2*drAbsInv
!    Gravity_bhEwaldAcc(1) = x*drAbsInv3 + grv_bhTreeEwald(4,i,j,k) + &
!    & x*grv_bhTreeEwald(7,i,j,k) + y*grv_bhTreeEwald(8,i,j,k) + z*grv_bhTreeEwald(9,i,j,k)
!    Gravity_bhEwaldAcc(2) = y*drAbsInv3 + grv_bhTreeEwald(5,i,j,k) + &
!    & x*grv_bhTreeEwald(8,i,j,k) + y*grv_bhTreeEwald(10,i,j,k) + z*grv_bhTreeEwald(11,i,j,k)
!    Gravity_bhEwaldAcc(3) = z*drAbsInv3 + grv_bhTreeEwald(6,i,j,k) + &
!    & x*grv_bhTreeEwald(9,i,j,k) + y*grv_bhTreeEwald(11,i,j,k) + z*grv_bhTreeEwald(12,i,j,k)

    Gravity_bhEwaldAcc(1) = grv_bhTreeEwald(4,i,j,k) + x*(drAbsInv3 + grv_bhTreeEwald(7,i,j,k)) + &
    & y*grv_bhTreeEwald(8,i,j,k) + z*grv_bhTreeEwald(9,i,j,k)
    Gravity_bhEwaldAcc(2) = grv_bhTreeEwald(5,i,j,k) +  &
    & x*grv_bhTreeEwald(8,i,j,k) + y*(drAbsInv3 + grv_bhTreeEwald(10,i,j,k)) + z*grv_bhTreeEwald(11,i,j,k)
    Gravity_bhEwaldAcc(3) = grv_bhTreeEwald(6,i,j,k) + &
    & x*grv_bhTreeEwald(9,i,j,k) + y*grv_bhTreeEwald(11,i,j,k) + z*(drAbsInv3 + grv_bhTreeEwald(12,i,j,k))


  endif


  return
end function
