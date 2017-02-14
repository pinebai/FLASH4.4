!!****if* source/Grid/GridSolvers/Multipole/gr_zoneMoments
!!
!! NAME
!!
!!  gr_zoneMoments
!!
!!
!! SYNOPSIS
!!    call gr_zoneMoments( real(IN) :: xprime,
!!                         real(IN) :: yprime,
!!                         real(IN) :: zprime,
!!                         real(IN) :: zonemass)
!!
!! DESCRIPTION
!!
!!    A routine to calculate the contribution of a cell to the
!!    moments of a mass distribution.
!!
!!    On exit, the current cell's contribution will have been added
!!    to Moment().  The associated Legendre polynomial values are
!!    calculated using a recurrence relation given in Press et al.
!!    (2d ed) pp 247-8.  The loops over l and m are arranged so as to
!!    most efficiently utilize this relation (ie. first loop over m,
!!    then over l).  The cosines and sines of multiples of the
!!    azimuthal angle are computed using the recurrence relations for
!!    the cosine and sine functions (see Press et al. p173).
!!
!! ARGUMENTS
!!
!!     xprime  -   X Coordinates of the current zone, relative to
!!                                 the center of mass
!!     yprime  -   Y Coordinates of the current zone, relative to
!!                                 the center of mass
!!     zprime  -   Z Coordinates of the current zone, relative to
!!                                 the center of mass
!!     zonemass  -    The mass in the current zone
!!
!!***



subroutine gr_zoneMoments (xprime, yprime, zprime, zonemass)

  use Grid_data, only : gr_meshMe
  use Logfile_interface, ONLY : Logfile_stamp

  use gr_mpoleData, ONLY : G_3DCARTESIAN,G_1DSPHERICAL,G_2DCYLINDRICAL,G_3DAXISYMMETRIC,&
                         Outer, Inner, Even, Odd, mpole_geometry, &
                         legk1, legk2, Moment, qmax, rpower, dsinv, rprinv,&
                         sintable, mpole_mmax, mpole_lmax, costable

  implicit none



  real, intent(in) :: xprime, yprime, zprime, zonemass


  !               "Local" variables

  !                 rprime        The length of the current position vector
  !                 lprime        The projection of rprime into the xy-plane
  !                 qprime        The index q (for Moment()) at which the current
  !                                 cell falls
  !                 mtemp1-7      Temporary variables
  !                 l,m,h         Temporary indices
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

  integer :: l, m, h, mp1, mp2, mm1, lmaxm1, lm1
  integer :: qprime, qmin
  real    :: rprime, lprime
  real    :: mtemp1, mtemp2, mtemp3, mtemp4, mtemp5, mtemp6, mtemp7
  real    :: trigalpha, trigbeta, costheta
  real    :: Legnd1, Legnd2, Legndr
  real    :: rscaled
  character(len=124) :: str_buffer

  !===============================================================================

  !                       Calculate the length of the position vector r'
  !                       and of its projection into the xy-plane (l').

  lprime = xprime**2 + yprime**2
  rprime = lprime + zprime**2
  lprime = sqrt(lprime)
  rprime = sqrt(rprime)

  call gr_mpoleEvalBin (rprime, rscaled)

  qprime = int(rscaled) + 1
  qmin   = 1



  if (qprime .gt. qmax) then
     write (str_buffer,*) 'gr_zoneMoments: desired q = ', qprime, & 
          ' is larger than qmax = ', qmax
     call Logfile_stamp( str_buffer, 'warning')
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

  Moment(qprime,Even,Inner,0,0) = Moment(qprime,Even,Inner,0,0) + rpower(0)
  Moment(qprime,Even,Outer,0,0) = Moment(qprime,Even,Outer,0,0) + rprinv(0)

  ! For 1D spherically symmetric problems only l = m = 0 contributes.

  if (mpole_geometry == G_1DSPHERICAL) return

  ! For 2D/3D problems, continue with l > 0.

  !                               Do (1,0).

  if (mpole_lmax >= 1) then
     Legnd2 = costheta * Legnd1
     mtemp3 = rpower(1) * Legnd2
     mtemp4 = rprinv(1) * Legnd2

     Moment(qprime,Even,Inner,1,0) = Moment(qprime,Even,Inner,1,0) + mtemp3
     Moment(qprime,Even,Outer,1,0) = Moment(qprime,Even,Outer,1,0) + mtemp4
  endif

  !                               Do (2,0) ... (mpole_lmax,0).

  do l = 2, mpole_lmax
     Legndr = costheta * Legk1(l,0) * Legnd2 - Legk2(l,0) * Legnd1
     Legnd1 = Legnd2
     Legnd2 = Legndr
     mtemp3 = rpower(l) * Legndr
     mtemp4 = rprinv(l) * Legndr

     Moment(qprime,Even,Inner,l,0) = Moment(qprime,Even,Inner,l,0) + mtemp3
     Moment(qprime,Even,Outer,l,0) = Moment(qprime,Even,Outer,l,0) + mtemp4
  enddo

  ! For 2D axisymmetric problems only m = 0 contributes.

  if (mpole_geometry == G_2DCYLINDRICAL) return

  ! also for 3D axisymmetric

  if (mpole_geometry == G_3DAXISYMMETRIC) return

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

     Moment(qprime,Even,Inner,m,m) = Moment(qprime,Even,Inner,m,m) + mtemp3
     Moment(qprime,Odd,Inner,m,m) = Moment(qprime,Odd,Inner,m,m) + mtemp4
     Moment(qprime,Even,Outer,m,m) = Moment(qprime,Even,Outer,m,m) + mtemp5
     Moment(qprime,Odd,Outer,m,m) = Moment(qprime,Odd,Outer,m,m) + mtemp6


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

     Moment(qprime,Even,Inner,mp1,m) = Moment(qprime,Even,Inner,mp1,m) + mtemp3
     Moment(qprime,Odd,Inner,mp1,m)  = Moment(qprime,Odd,Inner,mp1,m) + mtemp4
     Moment(qprime,Even,Outer,mp1,m) = Moment(qprime,Even,Outer,mp1,m) + mtemp5
     Moment(qprime,Odd,Outer,mp1,m)  = Moment(qprime,Odd,Outer,mp1,m) + mtemp6

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

        Moment(qprime,Even,Inner,l,m) = Moment(qprime,Even,Inner,l,m) + mtemp3
        Moment(qprime,Odd,Inner,l,m)  = Moment(qprime,Odd,Inner,l,m) + mtemp4
        Moment(qprime,Even,Outer,l,m) = Moment(qprime,Even,Outer,l,m) + mtemp5
        Moment(qprime,Odd,Outer,l,m)  = Moment(qprime,Odd,Outer,l,m) + mtemp6
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

     Moment(qprime,Even,Inner,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,Even,Inner,mpole_lmax,mpole_lmax) + mtemp3

     Moment(qprime,Odd,Inner,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,Odd,Inner,mpole_lmax,mpole_lmax) + mtemp4

     Moment(qprime,Even,Outer,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,Even,Outer,mpole_lmax,mpole_lmax) + mtemp5

     Moment(qprime,Odd,Outer,mpole_lmax,mpole_lmax) = & 
          Moment(qprime,Odd,Outer,mpole_lmax,mpole_lmax) + mtemp6
  endif

  !===============================================================================

  return
end subroutine gr_zoneMoments
