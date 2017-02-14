!!****f* source/physics/Gravity/Gravity_bhEwaldAccV42
!!
!! NAME
!!
!!  Gravity_bhEwaldAccV42
!!
!!
!! SYNOPSIS
!!
!!   real field = Gravity_bhEwaldAccV42(
!!                           real(in)    :: x,
!!                           real(in)    :: y,
!!                           real(in)    :: z
!!        )
!!
!! DESCRIPTION
!!
!!   Interpolates in the Ewald field and returns its value for the point x,y,z.
!!   Obsolete version from Flash4.2.
!!
!! ARGUMENTS
!!
!!  x   : x-coordinate of the point where the Ewald field is determined
!!  y   : y-coordinate of the point where the Ewald field is determined
!!  z   : z-coordinate of the point where the Ewald field is determined
!!
!! RESULT
!!
!!  Value of the Ewald field in a given point.
!!
!! NOTES
!!
!!***

#include "constants.h"

function Gravity_bhEwaldAccV42(x, y, z)
  implicit none
  real,intent(in) :: x, y, z
  real :: Gravity_bhEwaldAccV42(IAXIS:KAXIS)

  Gravity_bhEwaldAccV42 = 0.0

  return
end

