!!****f* source/physics/Cosmology/Cosmology_getOldRedshift
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
!!  the current cosmological scale factor based on the equation:
!!
!!    z = 1+(1/s)
!!   
!!  where 'z' is the last cosmolgical redshift, and 's' is the last
!!  scaling factor.
!!  
!!
!!***

subroutine Cosmology_getOldRedshift(zOld)

  implicit none
  real, intent(OUT) :: zOld

  zOld = 0.0

  return
end subroutine Cosmology_getOldRedshift
