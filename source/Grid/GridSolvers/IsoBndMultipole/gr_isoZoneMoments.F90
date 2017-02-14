!!****if* source/Grid/GridSolvers/IsoBndMultipole/gr_isoZoneMoments
!!
!! NAME
!!
!!  gr_isoZoneMoments
!!
!! SYNOPSIS
!!
!!  call gr_isoZoneMoments(real(IN) :: xprime,
!!                         real(IN) :: yprime,
!!                         real(IN) :: zprime,
!!                         real(IN) :: zonemass)
!!
!! DESCRIPTION
!!  A routine to calculate the contribution of a cell to the
!!    moments of a mass distribution.
!!
!!               On exit, the current cell's contribution will have been added
!!               to Moment().  The associated Legendre polynomial values are
!!               calculated using a recurrence relation given in Press et al.
!!               (2d ed) pp 247-8.  The loops over l and m are arranged so as to
!!               most efficiently utilize this relation (ie. first loop over m,
!!               then over l).  The cosines and sines of multiples of the
!!               azimuthal angle are computed using the recurrence relations for
!!               the cosine and sine functions (see Press et al. p173).
!!
!! ARGUMENTS
!!
!!   xprime, yprime, zprime -     Coordinates of the current zone, relative to
!!                                 the center of mass
!!   zonemass  -    The mass in the current zone
!!
!! NOTES
!!   The regular multipole solver (Grid/GridSolvers/Multipole) can handle
!!   the 3D axisymmetric case with the use of a parameter mpole_3daxisymmetric.
!!  This unit CANNOT, nor does it check for incorrect usage.
!!***



subroutine gr_isoZoneMoments (xprime, yprime, zprime, zonemass)

  !====================================================================

  use gr_isoMpoledata, ONLY: Moment, MPOLE_EVEN, MPOLE_ODD, MPOLE_INNER, MPOLE_OUTER,           &
       &     mpole_geometry, rpower, rprinv, sintable, costable, Legk1, Legk2, mpole_lmax, &
       &     qmax, dsinv, mpole_mmax

  use Logfile_interface, ONLY: Logfile_stampMessage

  implicit none
#include "Flash.h"
#include "constants.h"

  real, intent(IN) :: xprime, yprime, zprime, zonemass

  !               "Local" variables

  !                 rprime        The length of the current position vector
  !                 lprime        The projection of rprime into the xy-plane
  !                 qprime        The index q (for Moment()) at which the current
  !                                 cell falls
  !                 qlower, qupper
  !                               Limits (with qmin and qmax) for integrations
  !                 mtemp1-7      Temporary variables
  !                 l,m,h,q       Temporary indices
  !                 costable(m)   Table containing the cosine of m * the current
  !                                 azimuthal angle; also sintable(m)
  !                 trigalpha,
  !                 trigbeta      Trig coefficients for cosine/sine tables
  !                 costheta      The cosine of the current polar angle
  !                 Legnd1        The (m,m)th Legendre polynomial, evaluated at
  !                                 costheta
  !                 Legnd2        The (m+1,m)th Legendre polynomial
  !                 Legndr        The (l,m)th Legendre polynomial, with l>m+1
  !                 rpower(l)     The current cell's zonemass * rprime^l
  !                 rprinv(l)     The current cell's zonemass * rprime^-(l+1)
  !                 Leg_fact(l,m) Factorial normalization coefficients for the
  !                                 associated Legendre function

  integer :: l, m, h, q, mp1, mp2, mm1, lmaxm1, lm1
  integer :: qprime, qlower, qupper, qmin
  real    :: rprime, lprime
  real    :: mtemp1, mtemp2, mtemp3, mtemp4, mtemp5, mtemp6, mtemp7
  real    :: trigalpha, trigbeta, costheta
  real    :: Legnd1, Legnd2, Legndr

  character(len=128) :: str_buffer
  !===============================================================================

  !                       Calculate the length of the position vector r'
  !                       and of its projection into the xy-plane (l').

  lprime = xprime**2 + yprime**2
  rprime = lprime + zprime**2
  lprime = sqrt(lprime)
  rprime = sqrt(rprime)
  qprime = int(dsinv * rprime) + 1
  qmin   = 1
  qlower = max(qprime, qmin)
  qupper = min(qprime, qmax)

  !        print *, 'qmin, qmax, qlower, qupper = ', qmin, qmax, qlower, qupper

  if (qprime .gt. qmax) then
     write (str_buffer,*) 'gr_isoZonemoments:  WARNING:  desired q = ', qprime, & 
          ' is larger than qmax = ', qmax
     call Logfile_stampMessage(str_buffer)
  endif

  !                       Compute powers of r' * zonemass.

  mtemp1 = 1. / lprime
  mtemp2 = 1. / rprime
  rpower(0) = zonemass
  rprinv(0) = zonemass * mtemp2
  do l = 1, mpole_lmax
     lm1 = l - 1
     rpower(l) = rpower(lm1) * rprime
     rprinv(l) = rprinv(lm1) * mtemp2
  enddo

  !                       Compute table of cosines and sines for
  !                       multiples of the azimuthal angle phi';
  !                       also compute cos(theta').

  costable(0) = 1.
  sintable(0) = 0.
  trigalpha = 1. - (xprime * mtemp1)
  trigbeta = yprime * mtemp1
  do m = 1, mpole_mmax
     mm1 = m - 1
     costable(m) = costable(mm1) - & 
          (trigalpha*costable(mm1) + trigbeta*sintable(mm1))
     sintable(m) = sintable(mm1) - & 
          (trigalpha*sintable(mm1) - trigbeta*costable(mm1))
  enddo
  costheta = zprime * mtemp2
  mtemp7 = sqrt(1. - (costheta*costheta))

  !                       Compute the contributions to the
  !                       even and odd moments.

  !                       Do (l,0) moments.

  !                               Do (0,0).

  Legnd1 = 1.

  Moment(qprime,MPOLE_EVEN,MPOLE_INNER,0,0) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,0,0) + rpower(0)
  Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,0,0) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,0,0) + rprinv(0)

  ! For 1D spherically symmetric problems only l = m = 0 contributes.

  if ((NDIM==1).and.(mpole_geometry == SPHERICAL)) return

  ! For 2D/3D problems, continue with l > 0.

  !                               Do (1,0).

  if (mpole_lmax >= 1) then
     Legnd2 = costheta * Legnd1
     mtemp3 = rpower(1) * Legnd2
     mtemp4 = rprinv(1) * Legnd2

     Moment(qprime,MPOLE_EVEN,MPOLE_INNER,1,0) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,1,0) + mtemp3
     Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,1,0) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,1,0) + mtemp4
  endif

  !                               Do (2,0) ... (mpole_lmax,0).

  do l = 2, mpole_lmax
     Legndr = costheta * Legk1(l,0) * Legnd2 - Legk2(l,0) * Legnd1
     Legnd1 = Legnd2
     Legnd2 = Legndr
     mtemp3 = rpower(l) * Legndr
     mtemp4 = rprinv(l) * Legndr

     Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,0) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,0) + mtemp3
     Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,0) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,0) + mtemp4
  enddo

  ! For 2D axisymmetric problems only m = 0 contributes.

  if ((NDIM==2).and.(mpole_geometry == CYLINDRICAL)) return

! NOTE -- the regular gravity case returns with NDIM==3 & axisymmetric geometry
! NOTE     vecause the regular case can handle axisymmetric geometry in 2d & 3d, this routine can
! NOTE     can only tackle 2d axisymmetric

  ! For 3D problems, continue with m > 0.

  !                       Do (l,m) moments for 0 < m < mpole_lmax.

  lmaxm1 = mpole_lmax - 1
  do m = 1, lmaxm1

     !                               Do (m,m).

     Legnd1 = 1.
     mtemp5 = 1.
     do h = 1, m
        Legnd1 = -Legnd1 * mtemp5 * mtemp7
        mtemp5 = mtemp5 + 2.
     enddo
     mtemp3 = Legnd1 * costable(m)
     mtemp5 = rprinv(m) * mtemp3
     mtemp3 = rpower(m) * mtemp3
     mtemp4 = Legnd1 * sintable(m)
     mtemp6 = rprinv(m) * mtemp4
     mtemp4 = rpower(m) * mtemp4

     Moment(qprime,MPOLE_EVEN,MPOLE_INNER,m,m) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,m,m) + mtemp3
     Moment(qprime,MPOLE_ODD,MPOLE_INNER,m,m) = Moment(qprime,MPOLE_ODD,MPOLE_INNER,m,m) + mtemp4
     Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,m,m) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,m,m) + mtemp5
     Moment(qprime,MPOLE_ODD,MPOLE_OUTER,m,m) = Moment(qprime,MPOLE_ODD,MPOLE_OUTER,m,m) + mtemp6


     !                               Do (m+1,m).

     mp1 = m + 1
     mp2 = m + 2
     Legnd2 = costheta * (2*m+1) * Legnd1
     mtemp3 = Legnd2 * costable(m)
     mtemp5 = rprinv(mp1) * mtemp3
     mtemp3 = rpower(mp1) * mtemp3
     mtemp4 = Legnd2 * sintable(m)
     mtemp6 = rprinv(mp1) * mtemp4
     mtemp4 = rpower(mp1) * mtemp4

     Moment(qprime,MPOLE_EVEN,MPOLE_INNER,mp1,m) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,mp1,m) + mtemp3
     Moment(qprime,MPOLE_ODD,MPOLE_INNER,mp1,m)  = Moment(qprime,MPOLE_ODD,MPOLE_INNER,mp1,m) + mtemp4
     Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,mp1,m) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,mp1,m) + mtemp5
     Moment(qprime,MPOLE_ODD,MPOLE_OUTER,mp1,m)  = Moment(qprime,MPOLE_ODD,MPOLE_OUTER,mp1,m) + mtemp6

     !                               Do (m+2,m) ... (mpole_lmax,m).

     do l = mp2, mpole_lmax
        Legndr = costheta * Legk1(l,m) * Legnd2 - Legk2(l,m) * Legnd1
        Legnd1 = Legnd2
        Legnd2 = Legndr
        mtemp3 = Legndr * costable(m)
        mtemp5 = rprinv(l) * mtemp3
        mtemp3 = rpower(l) * mtemp3
        mtemp4 = Legndr * sintable(m)
        mtemp6 = rprinv(l) * mtemp4
        mtemp4 = rpower(l) * mtemp4

        Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,m) = Moment(qprime,MPOLE_EVEN,MPOLE_INNER,l,m) + mtemp3
        Moment(qprime,MPOLE_ODD,MPOLE_INNER,l,m)  = Moment(qprime,MPOLE_ODD,MPOLE_INNER,l,m) + mtemp4
        Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,m) = Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,l,m) + mtemp5
        Moment(qprime,MPOLE_ODD,MPOLE_OUTER,l,m)  = Moment(qprime,MPOLE_ODD,MPOLE_OUTER,l,m) + mtemp6
     enddo

  enddo

  !                               Do (mpole_lmax,mpole_lmax).

  if (mpole_lmax >= 1) then
     Legnd1 = 1.
     mtemp5 = 1.
     do h = 1, mpole_lmax
        Legnd1 = -Legnd1 * mtemp5 * mtemp7
        mtemp5 = mtemp5 + 2.
     enddo
     mtemp3 = Legnd1 * costable(mpole_lmax)
     mtemp5 = rprinv(mpole_lmax) * mtemp3
     mtemp3 = rpower(mpole_lmax) * mtemp3
     mtemp4 = Legnd1 * sintable(mpole_lmax)
     mtemp6 = rprinv(mpole_lmax) * mtemp4
     mtemp4 = rpower(mpole_lmax) * mtemp4

     Moment(qprime,MPOLE_EVEN,MPOLE_INNER,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,MPOLE_EVEN,MPOLE_INNER,mpole_lmax,mpole_lmax) + mtemp3

     Moment(qprime,MPOLE_ODD,MPOLE_INNER,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,MPOLE_ODD,MPOLE_INNER,mpole_lmax,mpole_lmax) + mtemp4

     Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,MPOLE_EVEN,MPOLE_OUTER,mpole_lmax,mpole_lmax) + mtemp5

     Moment(qprime,MPOLE_ODD,MPOLE_OUTER,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,MPOLE_ODD,MPOLE_OUTER,mpole_lmax,mpole_lmax) + mtemp6
  endif

  !===============================================================================

  return
end subroutine gr_isoZonemoments
