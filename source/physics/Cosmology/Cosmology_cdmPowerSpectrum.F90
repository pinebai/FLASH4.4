!!****f* source/physics/Cosmology/Cosmology_cdmPowerSpectrum
!!
!! NAME
!!
!!  Cosmology_cdmPowerSpectrum
!!
!! SYNOPSIS
!!
!!  call Cosmology_cdmPowerSpectrum(real, intent(IN)  :: k,
!!                                  real, intent(IN)  :: anorm,
!!                                  real, intent(IN)  :: znorm,
!!                                  real, intent(IN)  :: n,
!!                                  real, intent(IN)  :: omega0,
!!                                  real, intent(IN)  :: h,
!!                                  real, intent(IN)  :: lambda0,
!!                                  real, intent(OUT)  :: powerspectrum)
!!
!! DESCRIPTION
!!  
!!  Computes the present-day cold dark matter (CDM) power spectrum as a 
!!  function of a wavenumber.
!!
!! ARGUMENTS
!!
!!   k : Wavenumber
!!
!!   anorm : Cosmological scale factor normalization
!!
!!   znorm : Redshift normalization
!!
!!   n : Spectral index
!!
!!   omega0 : Present mass density
!!
!!   h : Hubble constant
!!
!!   lambda0 : Density parameter due to cosmological constant
!!
!!   powerspectrum : The generated power spectrum
!!
!!
!!
!!***

subroutine Cosmology_cdmPowerSpectrum(k, Anorm, znorm, n, Omega0, h, Lambda0, powerSpectrum)
    implicit none
    
    real, intent(IN) :: k, Anorm, znorm, n, Omega0, h, Lambda0
    real, intent(OUT) :: powerSpectrum

    powerSpectrum = 0.0
    return

end subroutine Cosmology_cdmPowerSpectrum
