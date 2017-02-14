!!****if* source/physics/sourceTerms/Flame/FlameSpeed/turbulent/tfi/charlette/gammaInt/fl_fsTFIGamma
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

  use fl_fsTFIInterface, ONLY : fl_fsTFIFunc
  use fl_fsTFIData, ONLY : fl_fsTFIus, fl_fsTFIdedl, fl_fsTFIPrandtl, &
                           fl_fsTFIetol

  implicit none

!  external fl_fsTFIExtFunc

  real, intent(OUT) :: GammaFit
  real, intent(IN) :: up_over_s, de_over_dl

  real :: result, a, b, abserr, etol2, etolabs, testf, kc

  real, parameter :: small_f = 1.e-8
  real, parameter :: large_k = 1.e15

  ! necessary work space for integrator
  integer, parameter :: limit = 100
  integer, parameter :: lenw = 4*limit
  integer, dimension(limit) :: iwork
  real, dimension(lenw) :: work
  integer :: key, last, neval, ier

  ! need to copy these to the module for the function to be integrated
  fl_fsTFIus = up_over_s
  fl_fsTFIdedl = de_over_dl

  ! integration is for Gamma^2, so etol should be modified
  etol2 = 2.0*fl_fsTFIetol
  etolabs = 0.0

  ! set limits for integration, a lower bound, b upper bound
  key = 3
  a = 1.0

  ! b should be infinity, but it dies to zero quickly
  ! guess from Eq. 26 from Charlette et al. (2002a)
  kc = min( 8.0/PI * (18.0*1.5/55.0)**1.5 * up_over_s**3, de_over_dl )
  b = max( 1.e3*kc, 10.0*a )
  testf = fl_fsTFIFunc(b)
  do while (testf > small_f)
     b = b*10.0
     if (b > large_k) &
        call Driver_abortFlash("fl_fsTFIGamma: Upper limit of integrand reached")
     testf = fl_fsTFIFunc(b)
  enddo

!  call dqag(fl_fsTFIExtFunc,a,b,etolabs,etol2,key,result,abserr, &
!           neval,ier,limit,lenw,last,iwork,work)
  call dqag(fl_fsTFIFunc,a,b,etolabs,etol2,key,result,abserr, &
           neval,ier,limit,lenw,last,iwork,work)

  if (ier.gt.0)  &
     call Driver_abortFlash("fl_fsTFIGamma: Some problem with the integrator")

  GammaFit = sqrt(result)

  return
end subroutine fl_fsTFIGamma
