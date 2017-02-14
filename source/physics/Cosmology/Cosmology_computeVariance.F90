!!****f* source/physics/Cosmology/Cosmology_computeVariance
!!
!! NAME
!!
!!  Cosmology_computeVariance
!!
!! SYNOPSIS
!!
!!  Cosmology_computeVariance(real,intent(inout)  :: lambda,
!!                            real,intent(inout)  :: mass,
!!                            real,intent(inout)  :: delta0,
!!                            real,intent(inout)  :: ddelta0dm,
!!                            integer, intent(in)  :: n,
!!                            real,intent(inout)  :: f,
!!                            real,intent(inout)  :: anorm,
!!                            real,intent(inout)  :: znorm,
!!                            real,intent(inout)  :: npspc,
!!                            real,intent(inout)  :: omega0,
!!                            real,intent(inout)  :: h,
!!                            real,intent(inout)  :: lambda0)
!!
!! DESCRIPTION
!!
!!  Given a range of comoving length scales, a processed power
!!  spectrum, a normalization redshift, and a smoothing kernel,
!!  compute the linear variance at the present epoch (ie.,
!!  (dM/M)^2).  The factor f is multiplied by the mass scale in
!!  applying the smoothing kernel.
!!
!! ARGUMENTS
!!
!!   lambda : Array of length scales
!!
!!   mass : Array of masses corresponding to length scales
!!
!!   delta0 : Arrays of linear variances (dM/M)^2
!!
!!   ddelta0dm : Deerivatives of linear variances with respect to mass
!!
!!   n : Number of elements in arrays
!!
!!   f : Factor of smoothing kernel
!!
!!   anorm : Cosmological scale factor normalization
!!
!!   znorm : Cosmological redshift normalization
!!
!!   npspc : Spectral index
!!
!!   omega0 : Present mass density
!!
!!   h : Hubble constant
!!
!!   lambda0 : Density parameter due to cosmological constant
!!
!!
!! NOTES
!!
!!  Uses the functions CDMPowerSpectrum and TopHatFilter from 
!!  MaterLambdaKernel/CosmologicalFunctions.F90.
!!
!!***

subroutine Cosmology_computeVariance (lambda, Mass, Delta0, dDelta0dM, N, f,&
                            Anorm, znorm, npspc, Omega0, h, Lambda0)

  implicit none

  integer, intent(in)  :: N
  real,intent(inout)     :: lambda(N), Delta0(N), f
  real,intent(inout)     :: dDelta0dM(N), Mass(N)
  real,intent(inout)     :: Anorm, znorm, npspc, Omega0, h, Lambda0
  
  
  return

end subroutine Cosmology_computeVariance
