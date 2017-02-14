!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoZonePotential
!!
!! NAME
!!
!!  gr_isoZonePotential
!!
!!
!! SYNOPSIS
!!
!!  call gr_isoZonePotential(real(IN) :: xprime, 
!!                      real(IN) :: yprime, 
!!                      real(IN) :: zprime, 
!!                      real(OUT):: potential)
!!
!! DESCRIPTION
!!
!!  A routine to calculate the value of the potential at one location using 
!!  the moments of the interior and exterior mass distributions.
!!
!!
!! ARGUMENTS
!!
!!  xprime, yprime, zprime  --  Coordinates of the current zone, relative to the center 
!!                     of mass
!!
!!  potential  --      The (estimated) value of the potential in the current 
!!                     zone (output)
!!
!! USES
!!
!!  From mpole_data
!!
!!  mpole_lmax          Maximum moment number (runtime parameter)
!!
!!  Moment(q,i,j,l,m)   The moments of the mass distribution (partially 
!!                      summed).  q = dividing radius; i = mpole_EVEN/MPOLE_ODD; 
!!                      j = MPOLE_INNER/MPOLE_OUTER; l, m = multipole indices.
!!
!!  costable(m)         Table containing the cosine of m * the current
!!                      azimuthal angle; also sintable(m)
!!
!!  rpower(l)           The value of rprime^l
!!
!!  rprinv(l)           The current cell's density * rprime^-(l+1)
!!
!!
!! NOTES
!!
!!  On exit, the potential at the specified location (actually, the estimate 
!!  based on the multipole expansion) is contained in potential.
!!
!!  The associated Legendre polynomial values are calculated using a recurrence 
!!  relation given in Press et al. (2d ed) pp 247-8.  The loops over l and m 
!!  are arranged so as to most efficiently utilize this relation (ie. first 
!!  loop over m, then over l).  The cosines and sines of multiples of the 
!!  azimuthal angle are computed using the recurrence relations for the cosine 
!!  and sine functions (see Press et al. p173).
!!
!!  To obtain the zone-averaged potential, mpole_potential calls this routine 
!!  for several different points within a zone, then finds the average of the 
!!  returned values.
!!
!!***

subroutine gr_isoZonepotential (xprime, yprime, zprime, potential)

#include "constants.h"
#include "Flash.h"

  use gr_isoMpoleData, ONLY: rprinv, dsinv, qmax, rpower, sintable, costable, &
       Legk1, Legk2, Moment, &
       mpole_lmax, mpole_mmax, mpole_geometry, &
       MPOLE_OUTER, MPOLE_INNER, MPOLE_EVEN, MPOLE_ODD

  use Logfile_interface, ONLY:  Logfile_stampMessage

  implicit none


  !            Local variables:
  !
  !               rprime          The length of the current position vector
  !               qprime          The index q (for Moment()) at which the current
  !                                 cell falls
  !               lprime          The projection of rprime into the xy-plane
  !               mtemp1-7         Temporary variables
  !               l,m,h,q         Temporary indices
  !               trigalpha,
  !               trigbeta        Trig coefficients for sine/cosine calculations
  !               costheta        The cosine of the current polar angle
  !               Legnd1          The (m,m)th Legendre polynomial, evaluated at
  !                                 costheta
  !               Legnd2          The (m+1,m)th Legendre polynomial
  !               Legndr          The (l,m)th Legendre polynomial, with l>m+1


  real, intent(IN)   :: xprime, yprime, zprime
  real, intent(OUT)  :: potential

  integer :: l, m, lm1, mp1, mm1, h, q, qprime, qprmp1, qprmm1
  real    :: mtemp1, mtemp2, mtemp3, mtemp4, mtemp5, mtemp7
  real    :: trigalpha, trigbeta, costheta
  real    :: Legnd1, Legnd2, Legndr, lprime, rprime
  real    :: rscaled, f, f1
  integer :: myPE

  character(len=128) :: str_buffer

  !===============================================================================

  !                       Clear the potential variable.

  potential = 0.

  !                       Calculate the length of the position vector r'
  !                       and of its projection into the xy-plane (l').

  lprime = xprime**2 + yprime**2
  rprime = lprime + zprime**2

  lprime = sqrt(lprime)
  rprime = sqrt(rprime)

  rscaled= rprime * dsinv

  qprime = int(rscaled) + 1
  qprmp1 = qprime + 1
  qprmm1 = qprime - 1

  f      = rscaled - qprmm1
  f1     = 1. - f

  if (qprime .gt. qmax) then
     write (str_buffer,*) 'gr_isoZonepotential:  WARNING:  desired q = ', qprime, & 
          ' is larger than qmax = ', qmax
     call Logfile_stampMessage(str_buffer)
  endif

  !                       Compute powers of r'.

  mtemp1 = 1. / lprime
  mtemp2 = 1. / rprime

  rpower(0) = 1.
  rprinv(0) = mtemp2

  ! Compute some factors for higher-order moments.

  if (mpole_lmax > 1) then

     do l = 1, mpole_lmax
        lm1 = l - 1
        rpower(l) = rpower(lm1) * rprime
        rprinv(l) = rprinv(lm1) * mtemp2
     enddo

     ! Compute table of cosines and sines for multiples of the azimuthal angle 
     ! phi'; also compute cos(theta').

     costable(0) = 1.
     sintable(0) = 0.

     trigalpha = 1. - (xprime * mtemp1)
     trigbeta = yprime * mtemp1

     do m = 1, mpole_mmax
        mm1 = m - 1
        costable(m) = costable(mm1) - (trigalpha*costable(mm1) + & 
             trigbeta*sintable(mm1))
        sintable(m) = sintable(mm1) - (trigalpha*sintable(mm1) - & 
             trigbeta*costable(mm1))
     enddo

     costheta = zprime * mtemp2
     mtemp7 = sqrt((1.-costheta)*(1.+costheta))

  endif

  !                       Compute the contributions of the even and odd moments
  !                       to the potential.

  !                       Do (l,0) moments for l = 0...mpole_lmax.

  !                               Do (0,0).

  Legnd1 = 1.
  potential = potential + & 
       (f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,0,0) + & 
       f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,0,0))*rprinv(0) + & 
       (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,0,0) + & 
       f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,0,0))*rpower(0)

  ! For 1D spherically symmetric problems only l = m = 0 contributes.

  if ((NDIM==1).and.(mpole_geometry == SPHERICAL)) return

  ! For 2D/3D problems, continue with l > 0.

  !                               Do (1,0).

  if (mpole_lmax >= 1) then
     Legnd2 = costheta * Legnd1
     potential = potential + & 
          (f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,1,0) + & 
          f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,1,0))*rprinv(1)*Legnd2 + & 
          (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,1,0) + & 
          f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,1,0))*rpower(1)*Legnd2
  endif

  !                               Do (2,0) ... (mpole_lmax,0).

  do l = 2, mpole_lmax
     Legndr = costheta * Legk1(l,0) * Legnd2 - Legk2(l,0) * Legnd1
     Legnd1 = Legnd2
     Legnd2 = Legndr
     potential = potential + & 
          (f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,l,0) + & 
          f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,0))*rprinv(l)*Legndr + & 
          (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,0) + & 
          f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,l,0))*rpower(l)*Legndr
  enddo

  ! For 2D axisymmetric problems only m = 0 contributes.

  if ((NDIM==2).and.(mpole_geometry == CYLINDRICAL)) return

!! NOTE regular Multipole gravity case returns with 3d axisymmetric as well, but this unit cannot
!!   handle that geometry

  ! For 3D problems, continue with m > 0.

  !                       Do (l,m) moments for m > 0.

  do m = 1, mpole_lmax

     mp1 = m + 1

     !                               Do (m,m).

     Legnd1 = 1.
     mtemp5 = 1.
     do h = 1, m
        Legnd1 = -Legnd1 * mtemp5 * mtemp7
        mtemp5 = mtemp5 + 2.
     enddo
     mtemp3 = rprinv(m) * Legnd1
     mtemp4 = rpower(m) * Legnd1
     potential = potential + & 
          ((f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,m,m) + & 
          f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,m,m))*mtemp3 + & 
          (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,m,m) + & 
          f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,m,m))*mtemp4)*costable(m) + & 
          ((f1*Moment(qprmm1,MPOLE_ODD,MPOLE_INNER,m,m) + & 
          f *Moment(qprime,MPOLE_ODD,MPOLE_INNER,m,m))*mtemp3 + & 
          (f1*Moment(qprime,MPOLE_ODD,MPOLE_OUTER,m,m) + & 
          f *Moment(qprmp1,MPOLE_ODD,MPOLE_OUTER,m,m))*mtemp4)*sintable(m)

     !                               Do (m+1,m).

     if (mp1 .le. mpole_lmax) then
        Legnd2 = costheta * (2*m+1) * Legnd1
        mtemp3 = rprinv(mp1) * Legnd2
        mtemp4 = rpower(mp1) * Legnd2
        potential = potential + & 
             ((f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,mp1,m) + & 
             f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,mp1,m))*mtemp3 + & 
             (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,mp1,m) + & 
             f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,mp1,m))*mtemp4)*costable(m) + & 
             ((f1*Moment(qprmm1,MPOLE_ODD,MPOLE_INNER,mp1,m) + & 
             f *Moment(qprime,MPOLE_ODD,MPOLE_INNER,mp1,m))*mtemp3 + & 
             (f1*Moment(qprime,MPOLE_ODD,MPOLE_OUTER,mp1,m) + & 
             f *Moment(qprmp1,MPOLE_ODD,MPOLE_OUTER,mp1,m))*mtemp4)*sintable(m)
     endif

     !                               Do (m+2,m) ... (mpole_lmax,m).

     do l = m+2, mpole_lmax
        Legndr = costheta * Legk1(l,m) * Legnd2 - Legk2(l,m) * Legnd1
        Legnd1 = Legnd2
        Legnd2 = Legndr
        mtemp3 = rprinv(l) * Legndr
        mtemp4 = rpower(l) * Legndr
        potential = potential + & 
             ((f1*Moment(qprmm1,MPOLE_EVEN,MPOLE_INNER,l,m) + & 
             f *Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,m))*mtemp3 + & 
             (f1*Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,m) + & 
             f *Moment(qprmp1,MPOLE_EVEN,MPOLE_OUTER,l,m))*mtemp4)*costable(m) + & 
             ((f1*Moment(qprmm1,MPOLE_ODD,MPOLE_INNER,l,m) + & 
             f *Moment(qprime,MPOLE_ODD,MPOLE_INNER,l,m))*mtemp3 + & 
             (f1*Moment(qprime,MPOLE_ODD,MPOLE_OUTER,l,m) + & 
             f *Moment(qprmp1,MPOLE_ODD,MPOLE_OUTER,l,m))*mtemp4)*sintable(m)
     enddo

  enddo


  return
end subroutine gr_isoZonepotential
