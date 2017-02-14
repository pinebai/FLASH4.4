!!****f* source/physics/Cosmology/Cosmology_computeDeltaCrit
!!
!! NAME
!!
!!  Cosmology_computeDeltaCrit
!!
!! SYNOPSIS
!!
!!  Cosmology_computeDeltaCrit(real, intent(in)  :: z,
!!                             real, intent(inout)  :: dcrit,
!!                             real, intent(inout)  :: dcritdz,
!!                             real, intent(inout)  :: d,
!!                             integer, intent(in)  :: n,
!!                             real, intent(in)  :: omega0,
!!                             real, intent(in)  :: h0,
!!                             real, intent(in)  :: lambda0,
!!                             real, intent(in)  :: omegatot)
!!
!! DESCRIPTION  
!!   
!!   Computes the (linear) overdensity at turnaround in the
!!   spherical collapse model at the given redshifts.  
!!   See Lacey and Cole 1993, MNRAS 262, 627 (appendix).
!!
!! ARGUMENTS
!!
!!   z : Array of redshifts
!!
!!   dcrit : Array of changes in the critical overdensity
!!
!!   dcritdz : Array of the critial overdensity derivative with respect to
!!             the cosmological redshift.
!!
!!   d : Array of linear perturbations at the given redshifts
!!
!!   n : Number of elements in arrays passed in 
!!
!!   omega0 : Present-day mass density
!!
!!   h0 : Hubble Constant
!!
!!   lambda0 : Cosmolofical constandt density 
!!
!!   omegatot : Total present-day mass density 
!!
!!
!!
!!***


subroutine Cosmology_computeDeltaCrit (z, dcrit, dcritdz, D, N, Omega0, H0, Lambda0, Omegatot)
  
  
  implicit none
  
  integer, intent(in) :: N
  real, intent(in)   :: z(N), Omega0, H0, Lambda0, Omegatot 
  real, intent(inout) :: D(N), dcrit(N), dcritdz(N)
  
  return

end subroutine Cosmology_computeDeltaCrit
