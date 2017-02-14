!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/charlette/gammaFit/fl_fsTFIGamma
!!
!! NAME
!!
!!  fl_fsTFIGamma
!!
!! SYNOPSIS
!!
!!  call fl_fsTFIGamma(real, intent(OUT)  :: gammafit,
!!                     real, intent(IN)  :: up_over_s,
!!                     real, intent(IN)  :: de_over_dl)
!!
!! DESCRIPTION
!!
!! Aaron Jackson 2010
!!
!! This subroutine calculates gamma
!!
!! ARGUMENTS
!!
!!   gammafit : 
!!
!!   up_over_s : 
!!
!!   de_over_dl : 
!!
!!
!!
!!***


#include "constants.h"

subroutine fl_fsTFIGamma(GammaFit, up_over_s, de_over_dl)

  use fl_fsTFIData, ONLY : fl_fsTFIPrandtl

  implicit none

  real, intent(OUT) :: GammaFit
  real, intent(IN) :: up_over_s, de_over_dl

  real, parameter :: Ck = 1.5
  real :: ReD, fu, fD, fRe, a, b, G

  ReD = up_over_s * de_over_dl / fl_fsTFIPrandtl

  ! calculate functions
  fu = 4.0*sqrt(27.0*Ck/110.0)*(18.0*Ck/55.0)*up_over_s**2
  fD = sqrt(27.0*Ck/110.0*PI**(4.0/3.0)*(de_over_dl**(4.0/3.0)-1.0))
  if (ReD > 0.0 ) then
    fRe = sqrt(9.0/55.0*exp(-1.5*Ck*PI**(4.0/3.0)/ReD)*ReD)
  else
    fRe = 0.0
  endif

  a = 0.60 + 0.20*exp(-0.1*up_over_s) - 0.20*exp(-0.01*de_over_dl)
  b = 1.4

  if (fRe > 0.0 ) then
    GammaFit = (((fu**(-a)+fD**(-a))**(-1.0/a))**(-b)+fRe**(-b))**(-1.0/b)
  else
    GammaFit = 0.0
  endif

  return
end subroutine fl_fsTFIGamma
