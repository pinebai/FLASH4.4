!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkDenseJakob
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
!!  routine networkDenseJakob sets up the dense aprox19 jacobian 
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

#include "Flash.h"

  
  use Burn_dataEOS
  use Burn_data
  use bn_dataAprox19

  implicit none

  ! arguments declaration
  integer, intent(IN) :: nlog, nphys
  real, intent(IN)    :: tt
  real, intent(OUT)   :: y(*), dfdy(nphys,nphys)

  ! Local variables
  integer          i,j
  real             yneut2,yprot2,xx,      &
       &                 den1,den2,r1,s1,t1,u1,v1,w1,x1,ralf1,ralf2,      &
       &                 r1f54,r2f54,r3f54,r4f54,r5f54,r6f54,r7f54,r8f54
  real             yy,den1a,den2a,r1f54a,r2f54a,                    &
       &                 r3f54a,r4f54a,r5f54a,r6f54a,r7f54a,r8f54a


  !..zero the jacobian
  real :: entropy, dst, dsd
  do j=1,NSPECIES
     do i=1,NSPECIES
        dfdy(i,j) = 0.0e0
     enddo
  enddo


  !..positive definite mass fractions
  do i=1,NSPECIES
     y(i) = min(1.0e0,max(y(i),1.0e-30))
  enddo



  !..branching ratios for dummy proton links
  !..and special combined rates for alpha photodisintegration and fe54
  yneut2 = y(ineut) * y(ineut)
  yprot2 = y(iprot) * y(iprot)
  r1     = ratdum(iralpa)/(ratdum(iralpa)+ratdum(iralpg))
  s1     = ratdum(irppa)/(ratdum(irppa)+ratdum(irppg))
  t1     = ratdum(irclpa)/(ratdum(irclpa)+ratdum(irclpg))
  u1     = ratdum(irkpa)/(ratdum(irkpa)+ratdum(irkpg))
  v1     = ratdum(irscpa)/(ratdum(irscpa)+ratdum(irscpg))
  w1     = ratdum(irvpa)/(ratdum(irvpa)+ratdum(irvpg))
  x1     = ratdum(irmnpa)/(ratdum(irmnpa)+ratdum(irmnpg))
  xx     = ratdum(irhegp) * ratdum(irdgn)                           &
       &        + y(ineut) * ratdum(irheng) * ratdum(irdgn)               &
       &        + y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)
  ralf1  = 0.0e0
  ralf2  = 0.0e0
  if (xx .ge. 1.0e-150) then
     ralf1 = ratdum(irhegn)*ratdum(irhegp)*ratdum(irdgn)/xx
     ralf2 = ratdum(irheng)*ratdum(irdpg)*ratdum(irhng)/xx
  end if


  !..don't consider 54fe photodisintegration if temperature is too low or there 
  !..is still free hydrogen around. (would like to use parameters here someday)
  xx  = ratdum(ircogp) + y(iprot)*(ratdum(ircopg) + ratdum(ircopa))
  if (xx .gt. 1.0e-99) then
     den1  = 1.0/xx
     den1a = -den1/xx * (ratdum(ircopg) + ratdum(ircopa))
  else
     den1  = 0.0e0
     den1a = 0.0e0
  end if
  if (btemp/1.0e9 .lt. 1.5 .or. y(ih1) .gt. 1.0e-5) then
     den1  = 0.0e0
     den1a = 0.0e0
  end if
  !
  yy  = ratdum(ir53gn) + y(ineut)*ratdum(ir53ng)
  if (yy .gt. 1.0e-99) then
     den2  = 1.0/yy
     den2a =  -den2/yy * ratdum(ir53ng)
  else
     den2  = 0.0e0
     den2a = 0.0e0
  end if
  if (btemp/1.0e9 .lt. 1.5 .or. y(ih1) .gt. 1.0e-5) then
     den2  = 0.0e0
     den2a = 0.0e0
  end if


  !       den1 = 0.0e0
  !       den2 = 0.0e0
  !       den1a = 0.0e0
  !       den2a = 0.0e0


  r1f54  = ratdum(ir54gn)*ratdum(ir53gn)*den2
  r1f54a = ratdum(ir54gn)*ratdum(ir53gn)*den2a
  r2f54  = ratdum(ir52ng)*ratdum(ir53ng)*den2
  r2f54a = ratdum(ir52ng)*ratdum(ir53ng)*den2a
  r3f54  = ratdum(irfepg)*ratdum(ircopg)*den1
  r3f54a = ratdum(irfepg)*ratdum(ircopg)*den1a
  r4f54  = ratdum(irnigp)*ratdum(ircogp)*den1
  r4f54a = ratdum(irnigp)*ratdum(ircogp)*den1a
  r5f54  = ratdum(irfepg)*ratdum(ircopa)*den1
  r5f54a = ratdum(irfepg)*ratdum(ircopa)*den1a
  r6f54  = ratdum(irfeap)*ratdum(ircogp)*den1
  r6f54a = ratdum(irfeap)*ratdum(ircogp)*den1a
  r7f54  = ratdum(irfeap)*ratdum(ircopg)*den1
  r7f54a = ratdum(irfeap)*ratdum(ircopg)*den1a
  r8f54  = ratdum(irnigp)*ratdum(ircopa)*den1
  r8f54a = ratdum(irnigp)*ratdum(ircopa)*den1a



  !..beta limit the he3+he4, 14n(pg), and 16o(pg) rates. assume steady state 
  !..abundances of 14o and 15o in beta limited cno cycle.
  if (y(ihe4) .gt. 0.0)                                             &
       &   ratdum(ir34)  = min(ratdum(ir34),0.896e0/y(ihe4))
  if (y(ih1)  .gt. 0.0)                                             &
       &   ratdum(irnpg) = min(ratdum(irnpg),5.68e-03/(y(ih1)*1.57e0))
  if (y(ih1)  .gt. 0.0)                                             &
       &   ratdum(iropg) = min(ratdum(iropg),0.0105e0/y(ih1))


  !..set up the jacobian
  !..hydrogen jacobian elements
  dfdy(ih1,ih1)   = -6.0e0 * y(ih1) * ratdum(irpp)                  &
       &                 - 2.0e0 * y(ic12) * ratdum(ircpg)                &
       &                 - 2.0e0 * y(in14) * ratdum(irnpg)                &
       &                 - 2.0e0 * y(io16) * ratdum(iropg)                &
       &                 - 3.0e0 * ratdum(irpen)
  dfdy(ih1,ihe3)  =  4.0e0 * y(ihe3) * ratdum(ir33)                 &
       &                 - y(ihe4) * ratdum(ir34) 
  dfdy(ih1,ihe4)  = -y(ihe3) * ratdum(ir34)
  dfdy(ih1,ic12)  = -2.0e0 * y(ih1) * ratdum(ircpg)
  dfdy(ih1,in14)  = -2.0e0 * y(ih1) * ratdum(irnpg)
  dfdy(ih1,io16)  = -2.0e0 * y(ih1) * ratdum(iropg)


  !..he3 jacobian elements
  dfdy(ihe3,ih1)  =  2.0e0 * ratdum(irpp)                           &
       &                 + ratdum(irpen)
  dfdy(ihe3,ihe3) = -4.0e0 * y(ihe3) * ratdum(ir33)                 &
       &                 - y(ihe4) * ratdum(ir34)
  dfdy(ihe3,ihe4) = -y(ihe3) * ratdum(ir34)


  !..he4 jacobian elements
  dfdy(ihe4,ih1)  =  y(in14) * ratdum(ifa) * ratdum(irnpg)          &
       &                 - y(io16) * ratdum(iropg)
  dfdy(ihe4,ihe3) =  2.0e0 * y(ihe3) * ratdum(ir33)                 &
       &                 - y(ihe4) * ratdum(ir34)
  dfdy(ihe4,ihe4) =  y(ihe3) * ratdum(ir34)                         &
       &                 - y(in14) * ratdum(irnag) * 1.5                  &
       &                 + r1 * y(img24) * ratdum(irmgap)                 &
       &                 + s1 * y(isi28) * ratdum(irsiap)                 &
       &                 + t1 * y(is32) * ratdum(irsap)                   &
       &                 + u1 * y(iar36) * ratdum(irarap)                 &
       &                 - y(io16) * ratdum(iroag)                        &
       &                 - y(ine20) * ratdum(irneag)                      &
       &                 - y(ic12) * ratdum(ircag)                        &
       &                 - y(img24) * (ratdum(irmgag)+ratdum(irmgap))     &
       &                 - y(isi28) * (ratdum(irsiag)+ratdum(irsiap))     &
       &                 - y(is32) * (ratdum(irsag)+ratdum(irsap))        &
       &                 - y(iar36) * (ratdum(irarag)+ratdum(irarap))     &
       &                 + v1 * y(ica40) * ratdum(ircaap)                 &
       &                 + w1 * y(iti44) * ratdum(irtiap)                 &
       &                 + x1 * y(icr48) * ratdum(ircrap)  
  dfdy(ihe4,ihe4) =  dfdy(ihe4,ihe4)                                &
       &                 - y(ife52) * y(iprot) * r7f54                    &
       &                 - y(ica40) * (ratdum(ircaag)+ratdum(ircaap))     &
       &                 - y(iti44) * (ratdum(irtiag)+ratdum(irtiap))     &
       &                 - y(icr48) * (ratdum(ircrag)+ratdum(ircrap))     &
       &                 - y(ife52) * (ratdum(irfeag)+r6f54)              &
       &                 - 9.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)       &
       &                 - ralf1
  dfdy(ihe4,ic12) =  2.0e0 * y(ic12) * ratdum(ir1212)               &
       &                 - y(ihe4) * ratdum(ircag)                        &
       &                 + 0.5e0 * y(io16) * ratdum(ir1216)               &
       &                 + 3.e0 * ratdum(irg3a) 
  dfdy(ihe4,in14) =  y(ih1) * ratdum(ifa) * ratdum(irnpg)           &
       &                 - y(ihe4) * ratdum(irnag) * 1.5
  dfdy(ihe4,io16) =  y(ih1) * ratdum(iropg)                         &
       &                 + 1.12e0 * y(io16) * ratdum(ir1616)              &
       &                 - y(ihe4) * ratdum(iroag)                        &
       &                 + 0.5e0 * y(ic12) * ratdum(ir1216)               &
       &                 + s1 * y(io16) * 0.68e0 * ratdum(ir1616)         &
       &                 + ratdum(iroga)
  dfdy(ihe4,ine20) = -y(ihe4) * ratdum(irneag)                      &
       &                  + ratdum(irnega)            
  dfdy(ihe4,img24) =  r1 * y(ihe4) * ratdum(irmgap)                 &
       &                  - y(ihe4) * (ratdum(irmgag)+ratdum(irmgap))     &
       &                 + ratdum(irmgga)
  dfdy(ihe4,isi28)=  s1 * y(ihe4) * ratdum(irsiap)                  &
       &                 - y(ihe4) * (ratdum(irsiag)+ratdum(irsiap))      &
       &                 + (ratdum(irsiga)+r1 * ratdum(irsigp)) 
  dfdy(ihe4,is32) =  t1 * y(ihe4) * ratdum(irsap)                   &
       &                 - y(ihe4) * (ratdum(irsag)+ratdum(irsap))        &
       &                  + (ratdum(irsga)+s1 * ratdum(irsgp))
  dfdy(ihe4,iar36) =  u1 * y(ihe4) * ratdum(irarap)                 &
       &                  - y(ihe4) * (ratdum(irarag)+ratdum(irarap))     &
       &                  + (ratdum(irarga) + t1*ratdum(irargp)) 
  dfdy(ihe4,ica40) =  v1 * y(ihe4) * ratdum(ircaap)                 &
       &                  - y(ihe4) * (ratdum(ircaag)+ratdum(ircaap))     &
       &                  + (ratdum(ircaga) + u1*ratdum(ircagp))
  dfdy(ihe4,iti44) =  w1 * y(ihe4) * ratdum(irtiap)                 &
       &                  - y(ihe4) * (ratdum(irtiag)+ratdum(irtiap))     &
       &                  + (ratdum(irtiga) + v1*ratdum(irtigp)) 
  dfdy(ihe4,icr48) =  x1 * y(ihe4) * ratdum(ircrap)                 &
       &                  - y(ihe4) * (ratdum(ircrag)+ratdum(ircrap))     &
       &                  + (ratdum(ircrga)+w1 * ratdum(ircrgp))
  dfdy(ihe4,ife52) = -y(ihe4) * y(iprot) * r7f54                    &
       &                  - y(ihe4) * (ratdum(irfeag)+r6f54)              &
       &                  + (ratdum(irfega)+x1 * ratdum(irfegp)) 
  dfdy(ihe4,ife54) =  yprot2 * r5f54
  dfdy(ihe4,ini56) =  (ratdum(irniga)+y(iprot) * r8f54)
  dfdy(ihe4,ineut) =  2.0e0 * y(ineut) * yprot2 * ralf2
  dfdy(ihe4,iprot) = -y(ife52) * y(ihe4)*r7f54 + y(ini56)*r8f54     &
       &                  + 2.0e0 * yneut2 * y(iprot) * ralf2             &
       &                  + 2.0e0 * y(ife54) * y(iprot) * r5f54           &
       &                  + y(ife54) * yprot2 * r5f54a                    &
       &                  - y(ihe4) * y(ife52) * r6f54a                   &
       &                  - y(ife52) * y(ihe4) * y(iprot) * r7f54a        &
       &                  + y(ini56) * y(iprot) * r8f54a


  !..c12 jacobian elements
  dfdy(ic12,ih1)    =  -y(ic12) * ratdum(ircpg)                     &
       &                    + y(in14) * ratdum(ifa) * ratdum(irnpg)
  dfdy(ic12,ihe4)   =  -y(ic12) * ratdum(ircag)                     &
       &                    + 3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
  dfdy(ic12,ic12)   =  -4.0e0 * y(ic12) * ratdum(ir1212)            &
       &                     - y(ihe4) * ratdum(ircag)                    &
       &                    - y(io16) * ratdum(ir1216)  - ratdum(irg3a)   &
       &                    - y(ih1) * ratdum(ircpg)          
  dfdy(ic12,in14)   =   y(ih1) * ratdum(ifa) * ratdum(irnpg)
  dfdy(ic12,io16)   =  -y(ic12) * ratdum(ir1216)   + ratdum(iroga)


  !..n14 jacobian elements
  dfdy(in14,ih1)    =  y(ic12) * ratdum(ircpg)                      &
       &                   - y(in14) * ratdum(irnpg)                      &
       &                   + y(io16) * ratdum(iropg)         
  dfdy(in14,ihe4)   = -y(in14) * ratdum(irnag)
  dfdy(in14,ic12)   =  y(ih1) * ratdum(ircpg)      
  dfdy(in14,in14)   = -y(ih1) * ratdum(irnpg)                       &
       &                   - y(ihe4) * ratdum(irnag)
  dfdy(in14,io16)   = y(ih1) * ratdum(iropg)


  !..16o jacobian elements
  dfdy(io16,ih1)    =  y(in14) * ratdum(ifg) * ratdum(irnpg)        &
       &                   - y(io16) * ratdum(iropg)
  dfdy(io16,ihe4)   = -y(io16) * ratdum(iroag)                      &
       &                   + y(ic12) * ratdum(ircag)
  dfdy(io16,ic12)   = -y(io16) * ratdum(ir1216)                     &
       &                   + y(ihe4) * ratdum(ircag)
  dfdy(io16,in14)   =  y(ih1) * ratdum(ifg) * ratdum(irnpg)      
  dfdy(io16,io16)   = -4.0e0 * y(io16) * ratdum(ir1616)             &
       &                   - y(ihe4) * ratdum(iroag)                      &
       &                   - y(ic12) * ratdum(ir1216)  - ratdum(iroga)    &
       &                   - y(ih1) * ratdum(iropg)
  dfdy(io16,ine20)  =  ratdum(irnega)


  !..20ne jacobian elements
  dfdy(ine20,ihe4)  = -y(ine20) * ratdum(irneag)                    &
       &                   + y(io16) * ratdum(iroag)                      &
       &                   + y(in14) * ratdum(irnag)
  dfdy(ine20,ic12)  =  2.0e0 * y(ic12) * ratdum(ir1212)      
  dfdy(ine20,in14)  =  y(ihe4) * ratdum(irnag)
  dfdy(ine20,io16)  =  y(ihe4) * ratdum(iroag)
  dfdy(ine20,ine20) = -y(ihe4) * ratdum(irneag) - ratdum(irnega)
  dfdy(ine20,img24) =  ratdum(irmgga)


  !..24mg jacobian elements
  dfdy(img24,ihe4)  = -y(img24)  *                                  &
       &                     (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1)) &
       &                   + y(ine20) * ratdum(irneag)
  dfdy(img24,ic12)  =  0.5e0 * y(io16) * ratdum(ir1216)
  dfdy(img24,io16)  =  0.5e0 * y(ic12) * ratdum(ir1216)
  dfdy(img24,ine20) =  y(ihe4) * ratdum(irneag)
  dfdy(img24,img24) = -y(ihe4)  *                                   &
       &                     (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1)) &
       &                   - ratdum(irmgga) 
  dfdy(img24,isi28) =  (ratdum(irsiga)+r1 * ratdum(irsigp))


  !..28si jacobian elements
  dfdy(isi28,ihe4)  =  y(img24)  *                                  &
       &                     (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1)) &
       &                   - y(isi28)  *                                  &
       &                     (ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1))
  dfdy(isi28,ic12)  =  0.5e0 * y(io16) * ratdum(ir1216) 
  dfdy(isi28,io16)  =  1.12e0 * y(io16) * ratdum(ir1616)            &
       &                   + 0.68e0 * y(io16) * s1 * ratdum(ir1616)       &
       &                   + 0.5e0 * y(ic12) * ratdum(ir1216) 
  dfdy(isi28,img24) =  y(ihe4)  *                                   &
       &                     (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1))
  dfdy(isi28,isi28) = -y(ihe4)  *                                   &
       &                     (ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1)) &
       &                   - (r1 * ratdum(irsigp)+ratdum(irsiga))
  dfdy(isi28,is32)  =  (ratdum(irsga)+s1 * ratdum(irsgp))


  !..32s jacobian elements
  dfdy(is32,ihe4)   =  y(isi28)  *                                  &
       &                     (ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1)) &
       &                   - y(is32)  *                                   &
       &                     (ratdum(irsag)+ratdum(irsap) * (1.0e0-t1))
  dfdy(is32,io16)   =  0.68e0 * y(io16)*ratdum(ir1616) * (1.0e0-s1) &
       &                   + 0.2e0 * y(io16) * ratdum(ir1616) 
  dfdy(is32,isi28)  =  y(ihe4)  *                                   &
       &                     (ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1))
  dfdy(is32,is32)   = -y(ihe4)  *                                   &
       &                     (ratdum(irsag)+ratdum(irsap) * (1.0e0-t1))   &
       &                   - (ratdum(irsga)+s1 * ratdum(irsgp))
  dfdy(is32,iar36)  =  (ratdum(irarga)+t1 * ratdum(irargp))


  !..36ar jacobian elements
  dfdy(iar36,ihe4)  =  y(is32)  *                                   &
       &                     (ratdum(irsag)+ratdum(irsap) * (1.0e0-t1))   &
       &                   - y(iar36)  *                                  &
       &                     (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1))
  dfdy(iar36,is32)  =  y(ihe4)  *                                   &
       &                     (ratdum(irsag)+ratdum(irsap) * (1.0e0-t1)) 
  dfdy(iar36,iar36) = -y(ihe4)  *                                   &
       &                     (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1)) &
       &                   - (ratdum(irarga) + t1 * ratdum(irargp))
  dfdy(iar36,ica40) =  (ratdum(ircaga) + ratdum(ircagp) * u1)


  !..40ca jacobian elements
  dfdy(ica40,ihe4)  =  y(iar36)  *                                  &
       &                     (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1)) &
       &                   - y(ica40)  *                                  &
       &                     (ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1))
  dfdy(ica40,iar36) =  y(ihe4)  *                                   &
       &                     (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1))
  dfdy(ica40,ica40) = -y(ihe4)  *                                   &
       &                     (ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1)) &
       &                   - (ratdum(ircaga) + ratdum(ircagp) * u1)
  dfdy(ica40,iti44) =  (ratdum(irtiga) + ratdum(irtigp) * v1)


  !..44ti jacobian elements
  dfdy(iti44,ihe4)  =  y(ica40)  *                                  &
       &                     (ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1)) &
       &                   - y(iti44)  *                                  &
       &                     (ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1))
  dfdy(iti44,ica40) =  y(ihe4)  *                                   &
       &                     (ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1))
  dfdy(iti44,iti44) = -y(ihe4)  *                                   &
       &                     (ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1)) &
       &                   - (ratdum(irtiga)+v1 * ratdum(irtigp))
  dfdy(iti44,icr48) =  (ratdum(ircrga)+w1 * ratdum(ircrgp))


  !..48cr jacobian elements
  dfdy(icr48,ihe4)  =  y(iti44)  *                                  &
       &                     (ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1)) &
       &                   - y(icr48)  *                                  &
       &                     (ratdum(ircrag)+ratdum(ircrap) * (1.0e0-x1))
  dfdy(icr48,iti44) =  y(ihe4)  *                                   &
       &                     (ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1))
  dfdy(icr48,icr48) = -y(ihe4)  *                                   &
       &                     (ratdum(ircrag)+ratdum(ircrap) * (1.0e0-x1)) &
       &                   - (ratdum(ircrga)+w1 * ratdum(ircrgp))
  dfdy(icr48,ife52) =  (ratdum(irfega)+x1 * ratdum(irfegp))


  !..52fe jacobian elements
  dfdy(ife52,ihe4)  =  y(icr48)  *                                  &
       &                     (ratdum(ircrag)+(1.0e0-x1) * ratdum(ircrap)) &
       &                   - y(ife52) * (ratdum(irfeag)+r6f54)            &
       &                   - y(ife52) * y(iprot) * r7f54  
  dfdy(ife52,icr48) =  y(ihe4)  *                                   &
       &                     (ratdum(ircrag)+(1.0e0-x1) * ratdum(ircrap))
  dfdy(ife52,ife52) = -(ratdum(irfega)+x1 * ratdum(irfegp))         &
       &                   - y(ihe4) * (ratdum(irfeag)+r6f54)             &
       &                   - y(ihe4) * y(iprot) * r7f54                   &
       &                   - yneut2 * r2f54
  dfdy(ife52,ife54) =  yprot2 * r5f54          + r1f54
  dfdy(ife52,ini56) =  ratdum(irniga)        + y(iprot)*r8f54
  dfdy(ife52,ineut) = -2.0e0 * y(ife52) * y(ineut) * r2f54          &
       &                     + y(ife54) * r1f54a                          &
       &                     - y(ife52) * yneut2 * r2f54a
  dfdy(ife52,iprot) =  2.0e0 * y(ife54) * y(iprot) * r5f54          &
       &                   - y(ife52)*y(ihe4)*r7f54  + y(ini56)*r8f54     &
       &                   + y(ife54) * yprot2 * r5f54a                   &
       &                   - y(ife52) * y(ihe4) * r6f54a                  &
       &                   - y(ife52) * y(ihe4) * y(iprot) * r7f54a       &
       &                   + y(ini56) * y(iprot) * r8f54a


  !..54fe jacobian elements
  dfdy(ife54,ihe4)  =  y(ife52) * r6f54  
  dfdy(ife54,ife52) =  yneut2 * r2f54    + y(ihe4) * r6f54           
  dfdy(ife54,ife54) = -r1f54           - yprot2 * (r3f54+r5f54) 
  dfdy(ife54,ini56) =  r4f54  + 56.0e0 * ratdum(irn56ec)/54.0e0
  dfdy(ife54,ineut) =  2.0e0 * y(ife52) * y(ineut) * r2f54          &
       &                     - y(ife54) * r1f54a                          &
       &                     + y(ife52) * yneut2 * r2f54a
  dfdy(ife54,iprot) = -2.0e0 * y(ife54) * y(iprot) * (r3f54+r5f54)  &
       &                   - y(ife54) * yprot2 * (r3f54a+r5f54a)          &
       &                    + y(ini56) * r4f54a                           &
       &                    + y(ife52) * y(ihe4) * r6f54a


  !..56ni jacobian elements
  dfdy(ini56,ihe4)  =  y(ife52) * (ratdum(irfeag)+y(iprot) * r7f54)
  dfdy(ini56,ife52) =  y(ihe4) * (ratdum(irfeag)+y(iprot) * r7f54)
  dfdy(ini56,ife54) =  yprot2 * r3f54  
  dfdy(ini56,ini56) = -(ratdum(irniga) + r4f54 + y(iprot) * r8f54)  &
       &                   - ratdum(irn56ec)
  dfdy(ini56,iprot) =  y(ife52) * y(ihe4) * r7f54  - y(ini56)*r8f54 &
       &                   + 2.0e0 * y(ife54) * y(iprot) * r3f54          &
       &                   + y(ife54) * yprot2 * r3f54a                   &
       &                   - y(ini56) * (r4f54a + y(iprot)*r8f54a)        &
       &                   + y(ife52) * y(ihe4) * y(iprot) * r7f54a


  !..neutron jacobian elements
  dfdy(ineut,ihe4)  =  2.0e0 * ralf1
  dfdy(ineut,ife52) = -2.0e0 * yneut2 * r2f54 
  dfdy(ineut,ife54) =  2.0e0 * r1f54        
  dfdy(ineut,ineut) =  -4.0e0 * y(ife52) * y(ineut) * r2f54         &
       &                    - 4.0e0 * yprot2 * y(ineut) * ralf2           &
       &                    + 2.0e0 * y(ife54) * r1f54a                   &
       &                    - 2.0e0 * y(ife52) * yneut2 * r2f54a          &
       &                    - ratdum(irnep)
  dfdy(ineut,iprot) = -4.0e0 * y(iprot) * yneut2 * ralf2            &
       &                     + ratdum(irpen)


  !..proton jacobian elements
  dfdy(iprot,ihe4)  =  2.0e0 * ralf1   + 2.0e0 * y(ife52) * r6f54
  dfdy(iprot,ife52) =  2.0e0 * y(ihe4) * r6f54
  dfdy(iprot,ife54) = -2.0e0 * yprot2 * (r3f54+r5f54)
  dfdy(iprot,ini56) =  2.0e0 * r4f54           
  dfdy(iprot,ineut) = -4.0e0 * yprot2 * y(ineut) * ralf2           &
       &                     + ratdum(irnep)
  dfdy(iprot,iprot) = -4.0e0 * y(iprot) * yneut2 * ralf2           &
       &                    - 4.0e0 * y(ife54) * y(iprot) * (r3f54+r5f54) &
       &                    - 2.0e0 * y(ife54) * yprot2 * (r3f54a+r5f54a) &
       &                    + 2.0e0 * y(ini56) * r4f54a                   &
       &                    + 2.0e0 * y(ife52) * y(ihe4) * r6f54a         &
       &                    - ratdum(irpen)

  return
end subroutine bn_networkDenseJakob





