!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_massToLength
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

  use Cosmology_data, ONLY : csm_omega, csm_hubble, csm_newton

  implicit none
  
  real, intent(IN) :: M
  real, intent(OUT) :: lambda
  real ::  M_arr(1), lambda_arr(1)

  M_arr(1) = M

  call MassToLengthConversion(M_arr, lambda_arr, 1, csm_omega, csm_hubble, csm_newton)

  lambda = lambda_arr(1)

  return
end subroutine Cosmology_massToLength



