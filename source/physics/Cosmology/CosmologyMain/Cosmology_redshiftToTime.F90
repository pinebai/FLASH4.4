!!***if* source/physics/Cosmology/CosmologyMain/Cosmology_redshiftToTime
!!
!! NAME
!!  Cosmology_redshiftToTime
!!
!! SYNOPSIS
!!  Cosmology_redshiftToTime(real,(IN)  :: z
!!                         real,(OUT) :: t)
!!
!! DESCRIPTION
!!
!!  Computes age of the universe corresponding to a given redshift for a given
!!  set of cosmological parameters.  The parameters involved are a part of 
!!  Cosmology_data and can be adjusted for a run in the flash.par file.
!!
!! ARGUMENTS
!!
!!  z -- A cosmological redshift value
!!  t -- The age of a universe at a given cosmological redshift
!!
!!***

subroutine Cosmology_redshiftToTime (z, t)
  
  use Cosmology_data, ONLY : csm_omega,csm_hubble,csm_lambda,csm_c
  
  implicit none

  real, intent(IN) :: z
  real, intent(OUT) :: t
  real :: Omegatot
  real :: z_arr(1), t_arr(1), dtdz(1)
  
  Omegatot = csm_omega + csm_lambda
  
  z_arr(1) = z
  
  call RedshiftToTimeConversion(z_arr, t_arr, dtdz, 1, csm_omega, &
       csm_hubble, csm_lambda, csm_c, Omegatot)

  t = t_arr(1)
  
  return
end subroutine Cosmology_redshiftToTime
