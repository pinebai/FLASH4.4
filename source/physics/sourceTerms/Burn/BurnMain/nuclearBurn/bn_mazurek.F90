!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bn_mazurek
!!
!! NAME
!!
!!   bn_mazurek
!!
!! SYNOPSIS
!!
!!   bn_mazurek(real(IN)  ::btemp,
!!              real(IN)  ::bden,
!!              real(IN)  ::y56,
!!              real(IN)  ::ye,
!!              real(OUT) ::rn56ec,
!!              real(OUT) ::sn56ec)
!!
!! DESCRIPTION
!!  routine mazurek computes ni56 electron capture rates
!!  this routine evaluates mazurel's 1973 fits for the ni56 electron
!!  capture rate rn56ec and neutrino loss rate sn56ec
!!
!! ARGUMENTS
!!
!!  btemp  - temperature -- routine returns 0.0 if (btemp .lt. 2.0e9)
!!  bden   - density -- routine returns 0.0 if (bden*ye .lt. 1.0e6)
!!  y56    - nickel56 molar abundance
!!  ye     - electron to baryon number, zbar/abar
!!  rn56ec - ni56 electron capture rate
!!  sn56ec - ni56 neutrino loss rate
!!
!!***


subroutine bn_mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)

  implicit none

  !!  declare arguments
  real, intent(IN) :: btemp, bden
  real, intent(IN) :: y56, ye
  real, intent(OUT):: rn56ec,sn56ec

  !!  declare local variables
  logical, save     ::   ifirst
  integer           :: j, k, jr, jd, ik, ij, jp, kp
  ! NOTE: LBR doesn't think the next line needs saving, but
  !  heck it's only a couple of scalars
  real, save        :: t9,r,rfm,rf0,rf1,rf2, &
       dfacm,dfac0,dfac1,dfac2,  &
       tfm,tf0,tf1,tf2,tfacm,tfac0,tfac1,tfac2
  real, save        ::   rnt(2),rne(2,7), &
       datn(2,6,7), tv(7),rv(6), &
       rfdm(4),rfd0(4),rfd1(4),rfd2(4),  &
       tfdm(5),tfd0(5),tfd1(5),tfd2(5)
  !!  initialize
  data  rv /6.e0, 7.e0, 8.e0, 9.e0, 10.e0, 11.e0/
  data  tv /2.e0, 4.e0, 6.e0, 8.e0, 10.e0, 12.e0, 14.e0/
  data ((datn(1,ik,ij),ik=1,6),ij=1,7) /  &
       -3.98e0, -2.84e0, -1.41e0,  0.20e0,  1.89e0,  3.63e0,  &
       -3.45e0, -2.62e0, -1.32e0,  0.22e0,  1.89e0,  3.63e0,  &
       -2.68e0, -2.30e0, -1.19e0,  0.27e0,  1.91e0,  3.62e0,  &
       -2.04e0, -1.87e0, -1.01e0,  0.34e0,  1.94e0,  3.62e0,  &
       -1.50e0, -1.41e0, -0.80e0,  0.45e0,  1.99e0,  3.60e0,  &
       -1.00e0, -0.95e0, -0.54e0,  0.60e0,  2.06e0,  3.58e0,  &
       -0.52e0, -0.49e0, -0.21e0,  0.79e0,  2.15e0,  3.55e0 /
  data ((datn(2,ik,ij),ik=1,6),ij=1,7) /  &
       -3.68e0, -2.45e0, -0.80e0,  1.12e0,  3.13e0,  5.19e0,  &
       -2.91e0, -2.05e0, -0.64e0,  1.16e0,  3.14e0,  5.18e0,  &
       -1.95e0, -1.57e0, -0.40e0,  1.24e0,  3.16e0,  5.18e0,  &
       -1.16e0, -0.99e0, -0.11e0,  1.37e0,  3.20e0,  5.18e0,  &
       -0.48e0, -0.40e0,  0.22e0,  1.54e0,  3.28e0,  5.16e0,  &
       0.14e0,  0.19e0,  0.61e0,  1.78e0,  3.38e0,  5.14e0,  &
       0.75e0,  0.78e0,  1.06e0,  2.07e0,  3.51e0,  5.11e0 /
  data  ifirst /.true./

  !!  first time; calculate the cubic interp parameters for ni56 electron capture
  if (ifirst) then
     do k=2,4
        rfdm(k)=1.e0/((rv(k-1)-rv(k))*(rv(k-1)-rv(k+1))*(rv(k-1)-rv(k+2)))
        rfd0(k)=1.e0/((rv(k)-rv(k-1))*(rv(k)-rv(k+1))*(rv(k)-rv(k+2)))
     enddo
     do k=2,4
        rfd1(k)=1.e0/((rv(k+1)-rv(k-1))*(rv(k+1)-rv(k))*(rv(k+1)-rv(k+2)))
        rfd2(k)=1.e0/((rv(k+2)-rv(k-1))*(rv(k+2)-rv(k))*(rv(k+2)-rv(k+1)))
     enddo
     do j=2,5
        tfdm(j)=1.e0/((tv(j-1)-tv(j))*(tv(j-1)-tv(j+1))*(tv(j-1)-tv(j+2)))
        tfd0(j)=1.e0/((tv(j)-tv(j-1))*(tv(j)-tv(j+1))*(tv(j)-tv(j+2)))
     enddo
     do j=2,5
        tfd1(j)=1.e0/((tv(j+1)-tv(j-1))*(tv(j+1)-tv(j))*(tv(j+1)-tv(j+2)))
        tfd2(j)=1.e0/((tv(j+2)-tv(j-1))*(tv(j+2)-tv(j))*(tv(j+2)-tv(j+1)))
     enddo
     ifirst = .false.
  end if

  !!  calculate ni56 electron capture and neutrino loss rates
  rn56ec = 0.e0
  sn56ec = 0.e0

  if ( (btemp .lt. 2.0e9) .or. (bden*ye .lt. 1.0e6)) return

  t9    = min(btemp,1.4e10) * 1.0e-9 
  r     = max(6.0e0,min(11.0e0,log10(bden*ye)))
  jp    = min(max(2,int(0.5e0*t9)),5)
  kp    = min(max(2,int(r)-5),4)
  rfm   = r - rv(kp-1)
  rf0   = r - rv(kp)
  rf1   = r - rv(kp+1)
  rf2   = r - rv(kp+2)
  dfacm = rf0*rf1*rf2*rfdm(kp)
  dfac0 = rfm*rf1*rf2*rfd0(kp)
  dfac1 = rfm*rf0*rf2*rfd1(kp)
  dfac2 = rfm*rf0*rf1*rfd2(kp)
  tfm   = t9 - tv(jp-1)
  tf0   = t9 - tv(jp)
  tf1   = t9 - tv(jp+1)
  tf2   = t9 - tv(jp+2)
  tfacm = tf0*tf1*tf2*tfdm(jp)
  tfac0 = tfm*tf1*tf2*tfd0(jp)
  tfac1 = tfm*tf0*tf2*tfd1(jp)
  tfac2 = tfm*tf0*tf1*tfd2(jp)

  !!  evaluate the spline fits
  do jr = 1,2
     do jd = jp-1,jp+2
        rne(jr,jd) = dfacm*datn(jr,kp-1,jd) &
             +dfac0*datn(jr,kp  ,jd) &
             +dfac1*datn(jr,kp+1,jd) &
             +dfac2*datn(jr,kp+2,jd)
     enddo
     rnt(jr) = tfacm*rne(jr,jp-1) &
          +tfac0*rne(jr,jp  ) &
          +tfac1*rne(jr,jp+1) &
          +tfac2*rne(jr,jp+2)
  enddo

  !!  set the output
  rn56ec = 10.0e0**rnt(1)
  sn56ec = 6.022548e+23 * 8.18683e-7 * y56 * 10.0e0**rnt(2)

  return

end subroutine bn_mazurek
