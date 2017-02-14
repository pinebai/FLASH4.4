!!****f* source/physics/Gravity/Gravity_bhEwaldPot
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

  implicit none
  real, intent(IN) :: x, y, z, drAbsInv
  real :: Gravity_bhEwaldPot

  Gravity_bhEwaldPot = 0.0

  return
end function
