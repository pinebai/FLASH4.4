!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_getOldRedshift
!!
!! NAME
!!
!!  Cosmology_getOldRedshift
!!
!! SYNOPSIS
!!
!!  Cosmology_getOldRedshift( real(OUT) :: zOld ) 
!!
!! ARGUMENTS
!! 
!!  zOld -- The simulation's last cosmological redshift
!!
!! DESCRIPTION
!!
!!  This routine returns the cosmological redshift of the last step based upon
!!  the last cosmological scale factor based on the equation:
!!
!!    z = 1+(1/s)
!!   
!!  where 'z' is the last cosmolgical redshift, and 's' is the last
!!  scaling factor.
!!  
!!
!!***

subroutine Cosmology_getOldRedshift(zOld)

  use Cosmology_data,  ONLY : csm_oldScaleFactor
  implicit none
  real, intent(OUT) :: zOld

  zOld = (1. / (csm_OldScaleFactor)) - 1.

  return
end subroutine Cosmology_getOldRedshift
