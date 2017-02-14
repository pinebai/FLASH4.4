!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_getParams
!!
!! NAME
!!
!!  Cosmology_getParams
!!
!! SYNOPSIS
!!
!!  Cosmology_getParams(real(OUT) :: hubble,
!!                      real(OUT) :: omega, 
!!                      real(OUT) :: baryon,
!!                      real(OUT) :: lambda) 
!!
!! ARGUMENTS
!! 
!!  hubble -- The simulation's value for the Hubble parameter.
!!  omega  -- The simulation's value for "Omega-Matter"
!!  baryon -- The simulation's value for "Omega-Baryon"
!!  lambda -- The simulation's value for "Omega-Lambda"
!!
!! DESCRIPTION
!!
!!  This routine returns the cosmological parameters used in the
!!  simulation.
!!
!!***

subroutine Cosmology_getParams(hubble,omega,baryon,lambda)

  use Cosmology_Data, ONLY: csm_hubble, csm_omega, csm_lambda, csm_baryon
  
  implicit none
  
  real, intent(OUT) :: hubble
  real, intent(OUT) :: omega
  real, intent(OUT) :: baryon
  real, intent(OUT) :: lambda

  hubble = csm_hubble
  omega = csm_omega
  baryon = csm_baryon
  lambda = csm_lambda

  return
end subroutine Cosmology_getParams
