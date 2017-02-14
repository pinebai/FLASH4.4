!*******************************************************************************

! Module:       CosmologicalFunctions

! Description:  This module provides a collection of functions that are useful
!               in initializing, performing, and analyzing cosmological
!               simulations.

! Provides:     MassToLength     Compute length scale given mass scale,
!                                obtaining cosmological parameters using the
!                                runtime parameter and physical constants
!                                databases.
!               MassToLengthConversion
!                                Compute length scales given mass scales for
!                                a given set of cosmological parameters.
!               CDMPowerSpectrum Compute present-day cold dark matter (CDM)
!                                power spectrum.
!               TopHatFilter     Tabulate the Fourier transform of the top-hat
!                                filter.
!               ComputeVariance  Compute the variance of a power spectrum.
!               RedshiftToTime   Compute universe age given redshift,
!                                obtaining cosmological parameters using the
!                                runtime parameter and physical constants
!                                databases.
!               RedshiftToTimeConversion
!                                Compute universe ages given redshifts for
!                                a given set of cosmological parameters.
!               ComputeDeltaCrit Compute linear overdensities at turnaround in
!                                the spherical collapse model.

! Planned:      TimeToRedshift
!               TimeToRedshiftConversion
!               LengthToMass
!               LengthToMassConversion
!               CriticalDensity
!               ***PowerSpectrum
!               ******Filter

! Local:        Integrate        Numerically integrate a tabulated function.
!               Integrand        Interpolate from a tabulated function to get
!                                an integrand value.
!               Convolve         Convolve two tabulated functions in k-space.
!               InitArray        Initialize an array with values spaced
!                                uniformly or logarithmically.


!===============================================================================

! Routine:      MassToLengthConversion

! Description:  Given a range of mass scales, compute the corresponding length
!               scales, ie. comoving diameters of spheres containing the given
!               amount of mass.


subroutine MassToLengthConversion (M, lambda, N, Omega0, H0, G)

  implicit none

  integer         :: N
  real            :: M(N), lambda(N)
  real            :: Omega0, G, H0
  
  real            :: K
  real, parameter :: onethd = 1./3.
  integer         :: i
  
  !-------------------------------------------------------------------------------
  
  K = ( 16.*G / (Omega0*H0**2) ) ** onethd
  do i = 1, N
     lambda(i) = K * M(i)**onethd
  enddo
  
  !-------------------------------------------------------------------------------
  
  return
end subroutine MassToLengthConversion

!===============================================================================

! Routine:      CDMPowerSpectrum

! Description:  Return the present-day CDM power spectrum at a given
!               wavenumber, given the normalization, the normalization
!               redshift, the primordial spectral index, Omega_0, and h.
!               The wavenumber and normalization must use length units of
!               Mpc, and the result is expressed in Mpc^3.


function CDMPowerSpectrum (k, Anorm, znorm, n, Omega0, h, Lambda0)
  
  implicit none
  
  real :: CDMPowerSpectrum
  real :: k, Anorm, znorm, n, Omega0, h, Lambda0
  
  real :: Omgh2, K1, K2, K3, K4, T, q
  
  !-------------------------------------------------------------------------------
  
  Omgh2 = 1. / (Omega0 * h**2 )
  K1    = Anorm                      ! should multiply by (D(0)/D(znorm))**2
  ! if znorm /= 0
  if (znorm /= 0.) &
       write (*,*) 'CDMPowerSpectrum:  WARNING:  znorm/=0 not yet implemented!'
  
  ! Kolb & Turner
  ! K2    = 1.7 * Omgh2
  ! K3    = 9.0 * Omgh2**1.5
  ! K4    = 1.0 * Omgh2**2
  ! CDMPowerSpectrum = K1 * k**(4+n) / (1. + K2*k + K3*k**1.5 + K4*k*k)**2
  
  ! Peebles
  ! K2 = 8. * Omgh2
  ! K3 = 4.7 * Omgh2**2
  ! CDMPowerSpectrum = K1 * k / (1. + K2*k + K3*k**2)**2
  
  ! BBKS eqn G3
  q = k * Omgh2
  T = log(1.+2.34*q) / (2.34*q) / &
       (1.+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**0.25
  CDMPowerSpectrum = K1 * k**(4+n) * T**2
  
  
  return
end function CDMPowerSpectrum

!===============================================================================

! Routine:      TopHatFilter

! Description:  Given a wavenumber and the characteristic length scale (cutoff
!               radius for the top hat function), return the Fourier transform
!               of the top hat filter.


function TopHatFilter (k, r)
  
  implicit none  
  real            :: TopHatFilter, k, r
  
  real            :: kr
  real, parameter :: pi = 3.141592654
  
  !-------------------------------------------------------------------------------
  
  kr = k * r
  TopHatFilter = 3. * (sin(kr) - kr*cos(kr)) / kr**3
  
  return
end function TopHatFilter

!===============================================================================

! Routine:      ComputeVariance

! Description:  Given a range of comoving length scales, a processed power
!               spectrum, a normalization redshift, and a smoothing kernel,
!               compute the linear variance at the present epoch (ie.,
!               (dM/M)^2).  The factor f is multiplied by the mass scale in
!               applying the smoothing kernel.


subroutine ComputeVariance (lambda, Mass, Delta0, dDelta0dM, N, f, PwrSpc, &
                            Filter, Anorm, znorm, npspc, Omega0, h, Lambda0)

  use Driver_interface, ONLY : Driver_abortFlash
  use ut_interpolationInterface

  implicit none

  integer  :: N
  real     :: lambda(N), Delta0(N), f, PwrSpc, Filter
  real     :: dDelta0dM(N), Mass(N)
  external    PwrSpc, Filter
  real     :: Anorm, znorm, npspc, Omega0, h, Lambda0
  
  real     :: pi, integral, kmin, kmax
  integer  :: Nint, i, j, m, mmax
  parameter   (Nint = 1000, kmin = 1.E-4, kmax = 1.E3)
  real     :: k(Nint), Pk(Nint), Wk(Nint)
  parameter   (pi = 3.141592654, mmax = 4)
  real     :: FuncP(mmax), FuncM(mmax), g, f3
  real     :: OrdP(mmax), OrdM(mmax)
  real     :: DerivF(mmax), DerivH(mmax), g3
  integer  :: ix,err
  
  !-------------------------------------------------------------------------------
  
  f3 = 0.5 * f**(1./3.)
  
  call InitArray (k, Nint, kmin, kmax, .true.)
  do j = 1, Nint
     Pk(j) = PwrSpc(k(j), Anorm, znorm, npspc, Omega0, h, Lambda0)
  enddo
  
  do j = 1, N
     
     do i = 1, Nint
        Wk(i) = Filter(k(i), lambda(j)*f3) ** 2      ! lambda is comoving diameter
     enddo
     
     call Convolve (k, Pk, Wk, Nint, integral)
     
     Delta0(j) = integral / (4.*pi*2.*pi**2)
     
     ! Use extrapolation to obtain the derivative.
     
     do m = 1, mmax
        g = 1. + 0.1/float(m**2)
        g3 = g**(1./3.)
        do i = 1, Nint
           Wk(i) = Filter(k(i), lambda(j)*f3*g3)**2
        enddo
        call Convolve (k, Pk, Wk, Nint, integral)
        FuncP(m) = integral / (4.*pi*2.*pi**2)
        OrdP(m)  = Mass(j)*f*g
        g = 1. - 0.1/float(m**2)
        g3 = g**(1./3.)
        do i = 1, Nint
           Wk(i) = Filter(k(i), lambda(j)*f3*g3)**2
        enddo
        call Convolve (k, Pk, Wk, Nint, integral)
        FuncM(m) = integral / (4.*pi*2.*pi**2)
        OrdM(m)  = Mass(j)*f*g
        DerivF(m) = (FuncP(m)-FuncM(m)) / (OrdP(m)-OrdM(m))
        DerivH(m) = OrdP(m)-OrdM(m)
     enddo
     
     call ut_fndpos(.FALSE.,DerivH,mmax,1,mmax,0.,ix,err)

     !call ut_polint (ix, DerivH, DerivF, mmax,  dDelta0dM(j), err)
     call Driver_abortFlash('ComputeVariance in CosmologicalFunctions: '//&
          'If this subroutine ever made sense, it certainly does not now!')
  enddo
  
  return
end subroutine ComputeVariance


!===============================================================================
! Routine:      RedshiftToTimeConversion

! Description:  Computes the age of the universe corresponding to each of a set
!               of given redshifts for a given set of cosmological parameters.


subroutine RedshiftToTimeConversion &
     (z, t, dtdz, N, Omega0, H0, Lambda0, c, Omegatot)
  
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer :: N
  real    :: z(N), t(N), dtdz(N), Omega0, H0, Lambda0, c, Omegatot
  
  integer :: i
  real    :: twothd, A, B, pi, eta, x
  parameter  (twothd = 2./3., pi = 3.14159265359)
  
  !-------------------------------------------------------------------------------
  
  if ((Omega0 .eq. 1.) .and. (Lambda0 .eq. 0.)) then      ! Flat universe
     
     do i = 1, N
        t(i)    = (twothd/H0) * (1.+z(i))**(-1.5)
        dtdz(i) = -1.5 * t(i) / (1.+z(i))
     enddo
     
  elseif ((Omega0 .lt. 1.) .and. (Lambda0 .eq. 0.)) then  ! Open universe
     
     B = Omega0 / (2.*(1.-Omega0)**1.5*H0)
     A = Omega0 / (2.*(1.-Omega0))
     do i = 1, N
        x    = 1./(A*(1.+z(i))) + 1.
        eta  = log( x + sqrt(x**2-1.) )                     ! eta = arccosh(x)
        ! log(x-sqrt(... for
        ! non-principal value
        t(i) = B * (sinh(eta) - eta)
        dtdz(i) = -A*B*(x-1.)**3/sqrt(x**2-1.)
     enddo
     
  elseif ((Omega0 .lt. 1.) .and. (Omegatot .eq. 1.)) then ! Lambda, flat
     
     do i = 1, N
        x    = sqrt((1.-Omega0)/Omega0) / (1.+z(i))**1.5
        t(i) = twothd/H0 / sqrt(1.-Omega0) * log( x + sqrt(x**2+1.) )
        ! arcsinh(x)
     enddo
     
     if (N .gt. 1) then
        do i = 2, N-1
           dtdz(i) = (t(i+1)-t(i-1)) / (z(i+1)-z(i-1))
        enddo
        dtdz(1) = (t(2)-t(1)) / (z(2)-z(1))
        dtdz(N) = (t(N)-t(N-1)) / (z(N)-z(N-1))
     else
        dtdz(1) = 0.
        write(*,*) 'WARNING:  RedshiftToTimeConversion:  dtdz not valid'
     endif
     
  else                                                    ! Unsupported
     
     call Driver_abortFlash("RedshiftToTimeConversion:  parameter values not yet implemented")

  endif
  
  return
end subroutine RedshiftToTimeConversion

!===============================================================================

! Routine:      ComputeDeltaCrit
! Description:  Computes the (linear) overdensity at turnaround in the
!               spherical collapse model at the given redshifts.  See Lacey and
!               Cole 1993, MNRAS 262, 627 (appendix).


subroutine ComputeDeltaCrit (z, dcrit, dcritdz, D, N, Omega0, H0, Lambda0, &
     Omegatot)
  
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  
  integer :: N
  real    :: z(N), dcrit(N), dcritdz(N), D(N), Omega0, H0, Lambda0, Omegatot
  
  integer :: i, Ny, j
  parameter  (Ny = 1000)
  real    :: y(Ny), yfunc(Ny)
  real    :: pi, twothd, d0crit, x, y1, y2, integral, D0
  parameter  (pi = 3.141592654, twothd = 2./3.)
  
  !-------------------------------------------------------------------------------
  
  if ((Omega0 .eq. 1.) .and. (Lambda0 .eq. 0.)) then       ! Flat universe
     ! Lacey & Cole 1993
     
     write (*,*) 'ComputeDeltaCrit:  flat universe, lambda=0'
     d0crit = 0.15 * (12.*pi)**twothd
     do i = 1, N
        D(i)       = 1. / (1.+z(i))
        dcrit(i)   = d0crit / D(i)
        dcritdz(i) = d0crit
     enddo
     
  elseif ((Omega0 .lt. 1.) .and. (Lambda0 .eq. 0.)) then   ! Open universe
     ! Henry 2000
     
     write (*,*) 'ComputeDeltaCrit:  open universe, lambda=0'
     x  = 2.*( (1./Omega0) - 1. )
     y1 = log(1.+x+sqrt(x**2+2*x))
     y2 = 1. + 6./x - 3.*sqrt(2.+x)*y1/x**1.5
     D0 = 1.25/((1./Omega0)-1.) * y2
     do i = 1, N
        x      = 2. * ( (1./Omega0) - 1. ) / (1.+z(i))
        y1     = log(1.+x+sqrt(x**2+2*x))
        y2     = 1. + 6./x - 3.*sqrt(2.+x)*y1/x**1.5
        d0crit = 1.5 * (1. + (2.*pi)**twothd / (sqrt(x**2+2*x) - y1)**twothd) * y2
        D(i)   = 1.25/((1./Omega0)-1.) * y2
        dcrit(i) = d0crit / D(i) * (D0 / (0.8*((1./Omega0)-1.)))
     enddo
     do i = 2, N-1
        dcritdz(i) = (dcrit(i+1)-dcrit(i-1)) / (z(i+1)-z(i-1))
     enddo
     dcritdz(1) = (dcrit(2)-dcrit(1)) / (z(2)-z(1))
     dcritdz(N) = (dcrit(N)-dcrit(N-1)) / (z(N)-z(N-1))
     
  elseif ((Omega0 .lt. 1.) .and.  (Omegatot .eq. 1.)) then ! Lambda, flat
     ! Henry 2000
     
     write (*,*) 'ComputeDeltaCrit:  flat universe, lambda/=0'
     call InitArray (y, Ny, 0., 1., .false.)
     y1 = ((1./Omega0)-1.)**(1./3.)
     do j = 1, Ny
        yfunc(j) = 1. / (1. + y1**3*y(j)**1.2)**1.5
     enddo
     call Integrate (y, yfunc, Ny, 0., 1., 1.E-10, 1.-1.E-10, integral)
     D0 = sqrt(1.+y1**3) * integral
     do i = 1, N
        x      = y1 / (1.+z(i))
        d0crit = 0.15*(12.*pi)**twothd * (1. - 0.0123*log10(1.+x**3))
        do j = 1, Ny
           yfunc(j) = 1. / (1. + x**3*y(j)**1.2)**1.5
        enddo
        call Integrate (y, yfunc, Ny, 0., 1., 1.E-10, 1.-1.E-10, integral)
        D(i) = (x/y1) * sqrt(1.+x**3) * integral
        dcrit(i) = d0crit / D(i) * D0
     enddo
     do i = 2, N-1
        dcritdz(i) = (dcrit(i+1)-dcrit(i-1)) / (z(i+1)-z(i-1))
     enddo
     dcritdz(1) = (dcrit(2)-dcrit(1)) / (z(2)-z(1))
     dcritdz(N) = (dcrit(N)-dcrit(N-1)) / (z(N)-z(N-1))
     
  else                                                     ! Unsupported
     
     call Driver_abortFlash ("ComputeDeltaCrit: unsupported parameter values")
     
  endif
  
  return
end subroutine ComputeDeltaCrit

!===============================================================================

! Subroutine:   Integrate

! Description:  Return the definite integral of a tabulated function.
!               Computes integral(x^p1*y(x)^p2, x=a..b) using the trapezoidal
!               rule.

! Parameters:   x       position samples
!               y       samples of function
!               N       number of samples
!               p1      power to raise position to
!               p2      power to raise function to
!               [a,b]   integration range
!               I       receives the value of the integral


subroutine Integrate (x, y, N, p1, p2, a, b, I)
  
  use Cosmology_data, ONLY: csm_meshMe
  use Driver_interface, ONLY : Driver_abortFlash

#include "constants.h"

  implicit none
  
  integer :: N
  real    :: x(N), y(N), a, b, I, p1, p2
  
  integer :: Ntab, j, nsteps, nstmax
  real    :: eps
  parameter  (Ntab = 1000, eps = 1.E-6, nstmax = 2**20)
  real    :: xtab(Ntab), ytab(Ntab), h, oldI, newsum
  
  common     /Integ/xtab, ytab
  real :: Integrand
  
  !-------------------------------------------------------------------------------
  
  if (N .gt. Ntab) call Driver_abortFlash ("Integrate:  N > Ntab!!!")
  
  do j = 1, N
     xtab(j) = x(j)
     ytab(j) = y(j)**p2 * x(j)**p1
  enddo
  do j = N+1, Ntab
     xtab(j) = x(N)
     ytab(j) = y(N)**p2 * x(N)**p1
  enddo
  
  nsteps = 1
  h      = b - a
  I      = 0.5*h*(Integrand(a) + Integrand(b))
  oldI   = 2.*I + 1.
  
  do while ( (nsteps < nstmax) .and. (abs(I-oldI) > eps*abs(I)) )
     oldI  = I
     h     = 0.5*h
     nsteps = nsteps*2
     newsum = 0.
     do j = 1, nsteps-1, 2
        newsum = newsum + Integrand(a+j*h)
     enddo
     I = 0.5*oldI + h*newsum
  enddo
  
  if (nsteps == nstmax) then
     if (csm_meshMe .EQ. MASTER_PE)print *, &
          'Integrate:  nsteps = ', nsteps, ' and error = ', &
          abs(I-OldI)/abs(I), ' > ', eps, '!'
  end if
  
  
  return
end subroutine Integrate



!===============================================================================

! Function:     Integrand

! Description:  Return the value of our tabulated integrand at some point.
!               Uses linear polynomial interpolation.


function Integrand (x)

  use Driver_interface, ONLY : Driver_abortFlash
  use ut_interpolationInterface

  implicit none
  
  real    :: Integrand, x
  real    :: y, err
  integer :: Ntab, j, k, m, ierr
  parameter  (Ntab = 1000, m = 2)
  real    :: xtab(Ntab), ytab(Ntab)
  common     /Integ/xtab, ytab
  
  !-------------------------------------------------------------------------------
  
  call ut_fndpos(.FALSE.,xtab,Ntab,1,Ntab,x,j,ierr)
  if (ierr .eq. -1) then
     call Driver_abortFlash('IERR from ut_fndpos says that something went wrong!')
  end if
  if ((j .eq. 0) .or. (j .eq. Ntab)) then
     write (*,*) 'Error - argument out of range in Integrand()'
     write (*,*) 'Argument = ', x
     write (*,*) 'Range    = [', xtab(1), '...', xtab(Ntab), ']'
     stop
  endif
  k = min( max(j-(m-1)/2, 1), Ntab+1-m )
  
  !!      y = ytab(j) + (ytab(j+1)-ytab(j))/(xtab(j+1)-xtab(j)) *
  !!     &              (x-xtab(j))
  call ut_polint (xtab(k), ytab(k), m, x, y, err)
  
  Integrand = y
  
  return
end function Integrand

!===============================================================================

! Routine:      Convolve
 
! Description:  Convolve two tabulated k-space functions, each of which depends
!               upon the magnitude of k only.  Assume the functions are zero for
!               |k| outside the tabulated range.  Integration is very crude
!               (trapezoidal rule).

 
subroutine Convolve (k, A, B, N, answer)
  
  implicit none
  
  integer :: N
  real    :: k(N), A(N), B(N), answer
  
  integer   :: Nmax, i, Nx
  parameter    (Nmax = 10000)
  real      :: x(Nmax), y(Nmax)
  
  !-------------------------------------------------------------------------------
  
  ! Copy the input data into local arrays so we can pass them to a function which
  ! interpolates from the table.
  
  if (N .gt. Nmax) then
     write (*,*) 'Convolve:  FATAL:  (N=',N,') > (Nmax=',Nmax,')'
     stop
  endif
  
  Nx = N
  do i = 1, Nx
     x(i) = k(i)
     y(i) = 4.*3.141592654 * x(i)**2 * A(i) * B(i)
  enddo
  
  ! Integrate using the trapezoidal rule.
  
  answer = 0.
  do i = 1, Nx-1
     answer = answer + 0.5*(y(i)+y(i+1))*(x(i+1)-x(i))
  enddo
  
 
  return
end subroutine Convolve

!===============================================================================

! Routine:      InitArray
! Description:  Initialize an array, optionally spacing values logarithmically.


subroutine InitArray (A, N, Amin, Amax, UseLog)

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  integer :: N
  real    :: A(N), Amin, Amax
  logical :: UseLog
  
  integer :: i
  real    :: dA, Alwr, Aupr
  
  !-------------------------------------------------------------------------------
  
  if (UseLog) then
     
     if ((Amin .le. 0.) .or. (Amax .le. 0.)) then
        call Driver_abortFlash("InitArray:  attempt to use log range with negative/zero limits")
     endif
     
     Alwr = log10(Amin)
     Aupr = log10(Amax)
     dA   = (Aupr - Alwr) / real(N-1)
     do i = 1, N
        A(i) = Amin * 10. ** (real(i-1)*dA)
     enddo
     
  else
     
     dA = (Amax - Amin) / real(N-1)
     do i = 1, N
        A(i) = Amin + real(i-1)*dA
     enddo
     
  endif
  
  return
end subroutine InitArray

