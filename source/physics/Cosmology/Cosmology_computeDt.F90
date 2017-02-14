!!****f* source/physics/Cosmology/Cosmology_computeDt
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
!!
!! ARGUMENTS
!!
!!  dt_cosmo -- variable to hold timestep constraint
!!
!!***


  
subroutine Cosmology_computeDt (dt_cosmo)

  implicit none
  
  real, INTENT(inout)    :: dt_cosmo

  return
end subroutine Cosmology_computeDt

