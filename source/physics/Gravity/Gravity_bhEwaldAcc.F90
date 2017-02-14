!!****f* source/physics/Gravity/Gravity_bhEwaldAcc
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

#include "constants.h"

function Gravity_bhEwaldAcc(x, y, z, drAbsInv)
  implicit none
  real, intent(IN) :: x, y, z, drAbsInv
  real :: Gravity_bhEwaldAcc(IAXIS:KAXIS)

  Gravity_bhEwaldAcc = 0.0

  return
end function
