!!****if* source/physics/Cosmology/CosmologyMain/csm_friedmannDeriv
!!
!! NAME
!!
!!  csm_friedmannDeriv
!!
!! SYNOPSIS
!!
!!  call csm_friedmannDeriv(real, intent(IN)  :: avar,
!!                          real, intent(INOUT)  :: dadtvar)
!!
!! DESCRIPTION
!!  
!!  This is the Friedmann Derivative equation that is used by
!!  csm_integrateFriedmann.  This will return the RHS of the 
!!  Friedmann equation.
!!
!! ARGUMENTS
!!
!!   avar :The initial scaling factor
!!
!!   dadtvar : the derivative of the scaling factor with respect to time
!!
!!
!!
!!***

subroutine csm_friedmannDeriv (avar, dadtvar)


  use Cosmology_data, ONLY : csm_hubble, csm_omega, csm_lambda, csm_curv
  implicit none
  
  real, intent(IN) :: avar
  real, intent(INOUT) :: dadtvar

  dadtvar = csm_hubble*sqrt(csm_omega/avar - csm_curv + csm_lambda*avar**2)

!===============================================================================

return
end subroutine csm_friedmannDeriv
