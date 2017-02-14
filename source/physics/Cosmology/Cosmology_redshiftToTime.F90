!!****f* source/physics/Cosmology/Cosmology_redshiftToTime
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

  implicit none
  
  real, intent(IN) :: z
  real, intent(OUT) :: t

  t = -1.0 !Set to a negative value.  Our universe has an age of <0 !?

end subroutine Cosmology_redshiftToTime
