!!****if* source/physics/Gravity/GravityMain/Poisson/BHTree/Gravity_bhEwaldPot
!!
!! NAME
!!
!!  Gravity_bhEwaldPot
!!
!! SYNOPSIS
!!
!!   real = Gravity_bhEwaldPot(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z,
!!                           real(in)    :: drAbsInv
!!        )
!!
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field for potential and returns its 
!!   value for the point x,y,z corrected to 1/r term
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
!!  Value of the normalised potential in a given point.
!!
!! NOTES
!! 
!!***

function Gravity_bhEwaldPot(x, y, z, drAbsInv)

  use Gravity_data, ONLY : grv_bhDxI, grv_bhTreeEwald, & 
      grv_bhPotConst, grv_margin, &
      grv_bhDLogI, grv_bhDLog, grv_bhLogfield, grv_bhLogfieldDer, &
      grv_bhEwald_periodicity, grv_bhLayerPotConst, grv_bhEwald_periodicity
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  real, intent(IN) :: x, y, z, drAbsInv

  integer :: i, j, k, jl
  real :: b, delta_jl, r2
  real :: Gravity_bhEwaldPot

  i = floor(0.5+x*grv_bhDxI)
  j = floor(0.5+y*grv_bhDxI)
  k = floor(0.5+z*grv_bhDxI)

! if the point is outside the Ewald field, analytical approximation is applied
  if ((max(i,j,k)).GT.grv_margin) then

    select case (grv_bhEwald_periodicity)
      case (1)
! evaluate the logarithm by interpolating in the grv_bhLogfield and grv_bhLogfieldDer tables
        r2 = y**2+z**2
        jl = floor(r2*grv_bhDLogI)
        delta_jl = r2 - jl*grv_bhDLog
!      Gravity_bhEwaldPot = (grv_bhLogfield(jl) + deltaR*grv_bhLogfieldDer(jl))*L1inv - grv_bhPotConst 
        Gravity_bhEwaldPot = grv_bhLogfield(jl) + delta_jl*grv_bhLogfieldDer(jl) + grv_bhPotConst 
      case (2)
        r2 = x**2+z**2
        jl = floor(r2*grv_bhDLogI)
        delta_jl = r2 - jl*grv_bhDLog
        Gravity_bhEwaldPot = grv_bhLogfield(jl) + delta_jl*grv_bhLogfieldDer(jl) + grv_bhPotConst
      case (4)
        r2 = x**2+y**2
        jl = floor(r2*grv_bhDLogI)
        delta_jl = r2 - jl*grv_bhDLog
        Gravity_bhEwaldPot = grv_bhLogfield(jl) + delta_jl*grv_bhLogfieldDer(jl) + grv_bhPotConst
      case (3)
        Gravity_bhEwaldPot = grv_bhLayerPotConst*z + grv_bhPotConst
      case (5)
        Gravity_bhEwaldPot = grv_bhLayerPotConst*y + grv_bhPotConst
      case (6)
        Gravity_bhEwaldPot = grv_bhLayerPotConst*x + grv_bhPotConst
      case (7)
        call Driver_abortFlash("This place must not be accessed")
    end select

  else
! interpolation and adding the 1/r term
    Gravity_bhEwaldPot =  drAbsInv + grv_bhTreeEwald(0,i,j,k) + x*grv_bhTreeEwald(1,i,j,k) + &
    &    y*grv_bhTreeEwald(2,i,j,k) + z*grv_bhTreeEwald(3,i,j,k)

  endif

  return
end function
