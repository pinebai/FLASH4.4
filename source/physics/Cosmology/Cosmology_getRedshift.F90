!!****f* source/physics/Cosmology/Cosmology_getRedshift
!!
!! NAME
!!
!!  Cosmology_getRedshift
!!
!! SYNOPSIS
!!
!!  Cosmology_getRedshift( real(OUT) :: z ) 
!!
!! ARGUMENTS
!! 
!!  z -- The simulation's current cosmological redshift
!!
!! DESCRIPTION
!!
!!  This routine returns the current cosmological redshift based upon
!!  the current cosmological scale factor based on the equation:
!!
!!    z = 1+(1/s)
!!   
!!  where 'z' is the current cosmolgical redshift, and 's' is the current
!!  scaling factor.
!!  
!!
!!***

subroutine Cosmology_getRedshift(z)

  implicit none
  real, intent(OUT) :: z

  z = 0.0

  return
end subroutine Cosmology_getRedshift
