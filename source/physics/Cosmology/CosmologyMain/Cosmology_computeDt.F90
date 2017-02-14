!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_computeDt
!!
!! NAME
!!
!!  Cosmology_computeDt
!!
!!
!! SYNOPSIS
!!
!!  Cosmology_computeDt(real,(INOUT) ::   dt_cosmo)
!!
!! DESCRIPTION
!!
!!  Computes the timestep limiter for the Cosmology Unit, based upon the 
!!  change in the scaling factor between timesteps, which is directly related
!!  to the change in the cosmological redshift.  
!!
!! ARGUMENTS
!!
!!  dt_cosmo --     variable to hold timestep constraint
!!
!!***


  
subroutine Cosmology_computeDt (dt_cosmo)

!===============================================================================

  use Cosmology_data, ONLY : csm_scaleFactor, csm_oldScaleFactor, &
       csm_maxScaleChange, csm_redshiftFinal
  use Driver_interface, ONLY : Driver_getDt, Driver_abortFlash

  implicit none
  
  real, INTENT(inout)    :: dt_cosmo
  
  real :: dadt, dt, afinal
  
!===============================================================================

  call Driver_getDt(dt)
 
  dadt = (csm_scaleFactor - csm_oldScaleFactor) / dt

  dt_cosmo = csm_maxScaleChange*csm_scaleFactor / dadt

! Try to avoid overshooting the final redshift
  afinal = 1./(1.+csm_redshiftFinal)

  if (csm_scaleFactor+2.*dadt*dt_cosmo > afinal) then
       dt_cosmo = (afinal - csm_scaleFactor) / (2.*dadt)
       dt_cosmo = MAX(dt_cosmo, TINY(1.0))
  endif


!===============================================================================

  if (dt_cosmo <= 0.0) call Driver_abortFlash("[Cosmology]: Computed dt is not positive! Aborting!")

  return
end subroutine Cosmology_computeDt

