!!****f* source/physics/Cosmology/Cosmology_massToLength
!!
!! NAME
!!
!!  Cosmology_massToLength
!!
!! SYNOPSIS
!!
!!  call Cosmology_massToLength(real, intent(IN)  :: m,
!!                              real, intent(OUT)  :: lambda)
!!
!! DESCRIPTION
!!
!!  Given a mass scale, compute the corresponding length scale, ie.the 
!!  comoving diameter of a sphere containing the given amount of mass.  
!!  Cosmological parameters are obtained from Cosmology_data
!!
!! ARGUMENTS
!!
!!   m : the mass scale
!!
!!   lambda : the corresponding length scale
!!
!!
!!
!!***

subroutine Cosmology_massToLength (M, lambda)

  implicit none
  
  real, intent(IN) :: M
  real, intent(OUT) :: lambda

  lambda = 0.0
  return
end subroutine Cosmology_massToLength
