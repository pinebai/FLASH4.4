!!****if* source/physics/Cosmology/CosmologyMain/Cosmology_cdmPowerSpectrum
!!
!! NAME
!!
!!  Cosmology_cdmPowerSpectrum
!!
!! SYNOPSIS
!!
!!  Cosmology_cdmPowerSpectrum(real(IN)  :: k, 
!!                             real(IN)  :: Anorm, 
!!                             real(IN)  :: znorm, 
!!                             real(IN)  :: n, 
!!                             real(IN)  :: Omega0, 
!!                             real(IN)  :: h, 
!!                             real(IN)  :: Lambda0,
!!                             real(OUT) :: powerSpectrum)
!!
!!
!! DESCRIPTION
!!  
!!  Computes the present-day cold dark matter (CDM) power spectrum as a 
!!  function of a given wavenumber.
!!
!! ARGUMENTS
!!
!!   k : Wavenumber
!!
!!   Anorm : Cosmological scale factor normalization
!!
!!   znorm : Redshift normalization
!!
!!   n : Spectral index
!!
!!   Omega0 : Present mass density
!!
!!   h : Hubble constant
!!
!!   Lambda0 : Density parameter due to cosmological constant
!!
!!   powerSpectrum : The generated power spectrum
!!
!!
!!
!!***


subroutine Cosmology_cdmPowerSpectrum(k, Anorm, znorm, n, Omega0, h, Lambda0, powerSpectrum)
    implicit none
    
    real, intent(IN) :: k, Anorm, znorm, n, Omega0, h, Lambda0
    real, intent(OUT) :: powerSpectrum
    
    real, external :: CDMPowerSpectrum

    powerSpectrum = CDMPowerSpectrum(k, Anorm, znorm, n, Omega0, h, Lambda0)

    return

end subroutine Cosmology_cdmPowerSpectrum
