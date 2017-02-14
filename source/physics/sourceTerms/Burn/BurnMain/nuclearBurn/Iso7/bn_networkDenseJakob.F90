!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_networkDenseJakob
!!
!! NAME
!!  
!!  bn_networkDenseJakob
!!
!! SYNOPSIS
!! 
!!  call bn_networkDenseJakob (real, intent(IN) :: tt,
!!                            real, intent(OUT) :: y(:),
!!                            real, intent(OUT) :: dfdy(nphys,nphys),
!!                          integer, intent(IN) :: nlog,
!!                          integer, intent(IN) :: nphys)
!!
!!  
!! DESCRIPTION
!!
!!  routine networkDenseJakob sets up the dense iso7 jacobian 
!!  
!!
!! ARGUMENTS
!!
!! tt :
!! y :
!! dfdy :
!! nlog :
!! nphys :
!!
!!***

subroutine bn_networkDenseJakob(tt,y,dfdy,nlog,nphys)   


  !      save
  use Burn_dataEOS, ONLY:  den, btemp
  use Burn_data
  use bn_dataIso7

  implicit none

#include "Flash.h"

  !
  !..this routine sets up the dense iso7 jacobian
  !
  !..declare arguments
  integer, intent(IN) :: nlog, nphys
  real, intent(IN)    :: tt
  real, intent(INOUT) ::  y(*)
  real, intent(OUT)   ::  dfdy(nphys,nphys)

  !..declare local variables
  integer    :: i,j
  real       :: t9,yeff_ca40,                &
       &           yeff_ti44,t9i,t932,t9i32,rsi2ni,rni2si,                &
       &           rsi2nida,rsi2nidsi,rni2sida


  !..zero the jacobian
  do j=1,nlog
     do i=1,nlog
        dfdy(i,j) = 0.0e0
     enddo
  enddo


  !..positive definite mass fractions
  do i=1,NSPECIES
     y(i) = min(1.0e0,max(y(i),1.0e-30))
  enddo


  !..some common factors and branching ratios
  t9    = btemp * 1.0e-9
  t9i   = 1.0e0/t9
  t932  = t9 * sqrt(t9)
  t9i32 = 1.0e0/t932


  !..rsi2ni is the rate for silicon to nickel
  !..rni2si is the rate for nickel to silicon
  rsi2ni    = 0.0e0
  rsi2nida  = 0.0e0
  rsi2nidsi = 0.0e0
  rni2si    = 0.0e0
  rni2sida  = 0.0e0

  !      if (t9 .gt. 20.0) then
  if (t9 .gt. 2.5 .and. y(ic12)+y(io16) .le. 4.0e-3) then
     yeff_ca40 = (t9i32**3) * exp(239.42*t9i-74.741)
     yeff_ti44 = (t932**3) * exp(-274.12*t9i+74.914)
     rsi2ni    = yeff_ca40*den**3*y(ihe4)**3*ratdum(ircaag)*y(isi28)
     rsi2nida  = 3.0e0 * rsi2ni/y(ihe4)
     rsi2nidsi = rsi2ni/y(isi28)
     rni2si = min(1.0e20,yeff_ti44*ratdum(irtiga)/(den**3*y(ihe4)**3))
     rni2sida  = -3.0e0 * rni2si/y(ihe4)
     if (rni2si .eq. 1.0e20) rni2sida = 0.0e0
  end if



  !..set up the jacobian
  !..4he jacobian elementss
  ! NOTE 2007/10/05 LBR.  First coefficient in front of y(ihe4) used to be -6.0, which 
  !   was wrong.  Problem noted by Chris Malone, confirmed by Frank Timmes.  Apparently
  !   doesn't happen in any other network.  
  dfdy(ihe4,ihe4)  = -9.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)      &
       &                  - y(ic12) * ratdum(ircag)                       &
       &                  - y(io16) * ratdum(iroag)                       &
       &                  - y(ine20) * ratdum(irneag)                     &
       &                  - y(img24) * ratdum(irmgag)                     &
       &                  - 7.0e0 * rsi2ni                                &
       &                  - 7.0e0 * rsi2nida * y(ihe4)                    &
       &                  + 7.0e0 * rni2sida * y(ini56)

  dfdy(ihe4,ic12)  =  3.0e0 * ratdum(irg3a)                         &
       &                  - y(ihe4) * ratdum(ircag)                       &
       &                  + 2.0e0 * y(ic12) * ratdum(ir1212)              &
       &                  + 0.5e0 * y(io16) * ratdum(ir1216)

  dfdy(ihe4,io16)  =  ratdum(iroga)                                 &
       &                  + 0.5e0 * y(ic12) * ratdum(ir1216)              &
       &                  + 2.0e0 * y(io16) * ratdum(ir1616)              &
       &                  - y(ihe4) * ratdum(iroag)

  dfdy(ihe4,ine20) =  ratdum(irnega)                                &
       &                  - y(ihe4) * ratdum(irneag) 

  dfdy(ihe4,img24) =  ratdum(irmgga)                                &
       &                  - y(ihe4) * ratdum(irmgag)

  dfdy(ihe4,isi28) =  ratdum(irsiga)                                &
       &                  - 7.0e0 * rsi2nidsi * y(ihe4)

  dfdy(ihe4,ini56) =  7.0e0 * rni2si



  !..12c jacobian elements 
  dfdy(ic12,ihe4) =  3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)       &
       &                  - y(ic12) * ratdum(ircag)

  dfdy(ic12,ic12) = -ratdum(irg3a)                                  &
       &                 - y(ihe4) * ratdum(ircag)                        &
       &                 - 4.0e0 * y(ic12) * ratdum(ir1212)               &
       &                 - y(io16) * ratdum(ir1216)

  dfdy(ic12,io16) =  ratdum(iroga)                                  &
       &                 - y(ic12) * ratdum(ir1216)


  !..16o jacobian elements 
  dfdy(io16,ihe4)  =  y(ic12) * ratdum(ircag)                       &
       &                  - y(io16) * ratdum(iroag) 

  dfdy(io16,ic12)  =  y(ihe4) * ratdum(ircag)                       &
       &                  - y(io16) * ratdum(ir1216) 

  dfdy(io16,io16)  = -ratdum(iroga)                                 &
       &                  - y(ic12) * ratdum(ir1216)                      &
       &                  - 4.0e0 * y(io16) * ratdum(ir1616)              &
       &                  - y(ihe4) * ratdum(iroag) 

  dfdy(io16,ine20) =  ratdum(irnega)



  !..20ne jacobian elements 
  dfdy(ine20,ihe4)  =  y(io16) * ratdum(iroag)                      &
       &                   - y(ine20) * ratdum(irneag)

  dfdy(ine20,ic12)  =  2.0e0 * y(ic12) * ratdum(ir1212)

  dfdy(ine20,io16)  =  y(ihe4) * ratdum(iroag)

  dfdy(ine20,ine20) = -ratdum(irnega)                               &
       &                   - y(ihe4) * ratdum(irneag)

  dfdy(ine20,img24) =  ratdum(irmgga)



  !..24mg jacobian elements 
  dfdy(img24,ihe4)  =  y(ine20) * ratdum(irneag)                    &
       &                   - y(img24) * ratdum(irmgag)

  dfdy(img24,ic12)  =  0.5e0 * y(io16) * ratdum(ir1216)

  dfdy(img24,io16)  =  0.5e0 * y(ic12) * ratdum(ir1216)

  dfdy(img24,ine20) =  y(ihe4) * ratdum(irneag)  

  dfdy(img24,img24) = -ratdum(irmgga)                               &
       &                   - y(ihe4) * ratdum(irmgag)

  dfdy(img24,isi28) =  ratdum(irsiga)



  !..28si jacobian elements 
  dfdy(isi28,ihe4)  =  y(img24) * ratdum(irmgag)                    &
       &                   - rsi2ni                                       &
       &                   - rsi2nida * y(ihe4)                           &
       &                   + rni2sida * y(ini56)

  dfdy(isi28,ic12)  = 0.5e0 * y(io16) * ratdum(ir1216) 

  dfdy(isi28,io16)  =  2.0e0 * y(io16) * ratdum(ir1616)             &
       &                   + 0.5e0 * y(ic12) * ratdum(ir1216) 

  dfdy(isi28,img24) =  y(ihe4) * ratdum(irmgag) 

  dfdy(isi28,isi28) = -ratdum(irsiga)                               &
       &                   - rsi2nidsi * y(ihe4)

  dfdy(isi28,ini56) =  rni2si 



  !..ni56 jacobian elements 
  dfdy(ini56,ihe4)  =   rsi2ni                                      &
       &                    + rsi2nida * y(ihe4)                          &
       &                    - rni2sida * y(ini56)

  dfdy(ini56,isi28) =  rsi2nidsi * y(ihe4)

  dfdy(ini56,ini56) = -rni2si

  return
end subroutine bn_networkDenseJakob

