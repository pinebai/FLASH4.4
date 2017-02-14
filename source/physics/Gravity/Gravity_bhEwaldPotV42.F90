!!****f* source/physics/Gravity/Gravity_bhEwaldPotV42
!!
!! NAME
!!
!!  Gravity_bhEwaldPotV42
!!
!!
!! SYNOPSIS
!!
!!   real field = Gravity_bhEwaldPotV42(
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



function Gravity_bhEwaldPotV42(x, y, z)

  implicit none
  real,intent(in) :: x, y, z
  real :: Gravity_bhEwaldPotV42

  Gravity_bhEwaldPotV42 = 0.0


  return
end

