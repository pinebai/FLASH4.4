!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkSparseJakob
!!
!! NAME
!!
!!  bn_networkSparseJakob
!!
!! SYNOPSIS
!! 
!!  call bn_networkSparseJakob(real, intent(IN)    :: tt,
!!                             real, intent(IN)    :: y(:),
!!                             real, intent(OUT)   :: dfdy(nphys,nphys),
!!                             integer, intent(IN) :: nzo,
!!                             integer, intent(IN) :: nDummy)
!!
!!  
!! DESCRIPTION
!!
!!  routine networkSparseJakob sets up the sparse aprox13 jacobian 
!!  input is tt (irrelevant here) and the abundances y(*). 
!!  output is the jacobian dfdy(nzo).
!!
!!  This routine is one of (a few) that can be called as an external
!!    'jakob'
!!
!! ARGUMENTS
!!
!!   tt -- not used in this subroutine
!!   y  -- abundances
!!   dfdy -- jacobian
!!   nzo -- size of jacobian
!!   nDummy -- dummy argument to make argument lists consistent across networks.
!!
!!***
subroutine bn_networkSparseJakob(tt,y,dfdy,nzo,nDummy)

#include "Flash.h"

  use Driver_interface, ONLY : Driver_abortFlash

  use Burn_dataEOS, ONLY:  btemp
  use Burn_data
  use bn_dataNetworkSize, ONLY:  neloc, nterms, eloc
  use bn_dataAprox19

  implicit none

  !
  !..this routine sets up the sparse aprox19 jacobian. 
  !..input is tt (irrelevant here) and the abundances y(*). 
  !..output is the jacobian dfdy(nzo).


  !..declare
  integer, intent(IN) :: nzo, nDummy
  real, intent(IN)    :: tt, y(*)
  real, intent(OUT)   :: dfdy(*)

  integer i,nt,iat
  real    a1,yneut2,yprot2,xx,den1,den2,            &
       &        r1,s1,t1,u1,v1,w1,x1,ralf1,ralf2,                         &
       &        r1f54,r2f54,r3f54,r4f54,r5f54,r6f54,r7f54,r8f54,          &
       &        yy,den1a,den2a,r1f54a,r2f54a,                             &
       &        r3f54a,r4f54a,r5f54a,r6f54a,r7f54a,r8f54a,dum1,dum2


  !..zero the jacobian
  real :: entropy, dst, dsd
  nt = 0
  do i=1,nzo
     dfdy(i) = 0.0e0
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
  !      den1a = 0.0e0
  !      den2a = 0.0e0


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



  !..hydrogen jacobian elements
  !..d(h1)/d(h1)
  a1 = -6.0e0 * y(ih1) * ratdum(irpp)                               &
       &    - 2.0e0 * y(ic12) * ratdum(ircpg)                             &
       &    - 2.0e0 * y(in14) * ratdum(irnpg)                             &
       &    - 2.0e0 * y(io16) * ratdum(iropg)                             &
       &    - 3.0e0 * ratdum(irpen)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(h1)/d(he3)
  a1  =  4.0e0 * y(ihe3) * ratdum(ir33) - y(ihe4) * ratdum(ir34) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(h1)/d(he4)
  a1  = -y(ihe3) * ratdum(ir34)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(h1)/d(c12)
  a1  = -2.0e0 * y(ih1) * ratdum(ircpg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(h1)/d(n14)
  a1  = -2.0e0 * y(ih1) * ratdum(irnpg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(h1)/d(o16)
  a1  = -2.0e0 * y(ih1) * ratdum(iropg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1




  !..he3 jacobian elements
  !..d(he3)/d(h1)
  a1  =  2.0e0 * ratdum(irpp)  + ratdum(irpen)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he3)/d(he3)
  a1  = -4.0e0 * y(ihe3) * ratdum(ir33) - y(ihe4) * ratdum(ir34)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he3)/d(he4)
  a1  = -y(ihe3) * ratdum(ir34)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1




  !..he4 jacobian elements
  !..d(he4)/d(h1)
  a1  =  y(in14)*ratdum(ifa) * ratdum(irnpg) - y(io16)*ratdum(iropg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(he3)
  a1  = 2.0e0 * y(ihe3) * ratdum(ir33) - y(ihe4) * ratdum(ir34)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(he4)
  a1 =  y(ihe3) * ratdum(ir34)                                      &
       &    - y(in14) * ratdum(irnag) * 1.5e0                             &
       &    + r1 * y(img24) * ratdum(irmgap)                              &
       &    + s1 * y(isi28) * ratdum(irsiap)                              &
       &    + t1 * y(is32) * ratdum(irsap)                                &
       &    + u1 * y(iar36) * ratdum(irarap)                              &
       &    - y(io16) * ratdum(iroag)                                     &
       &    - y(ine20) * ratdum(irneag)                                   &
       &    - y(ic12) * ratdum(ircag)                                     &
       &    - y(img24) * (ratdum(irmgag)+ratdum(irmgap))                  &
       &    - y(isi28) * (ratdum(irsiag)+ratdum(irsiap))                  &
       &    - y(is32) * (ratdum(irsag)+ratdum(irsap))                     &
       &    - y(iar36) * (ratdum(irarag)+ratdum(irarap))                  &
       &    + v1 * y(ica40) * ratdum(ircaap)                              &
       &    + w1 * y(iti44) * ratdum(irtiap)                              &
       &    + x1 * y(icr48) * ratdum(ircrap)  
  a1 =  a1                                                          &
       &    - y(ife52) * y(iprot) * r7f54                                 &
       &    - y(ica40) * (ratdum(ircaag)+ratdum(ircaap))                  &
       &    - y(iti44) * (ratdum(irtiag)+ratdum(irtiap))                  &
       &    - y(icr48) * (ratdum(ircrag)+ratdum(ircrap))                  &
       &    - y(ife52) * (ratdum(irfeag)+r6f54)                           &
       &    - 9.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)                    &
       &    - ralf1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(c12)
  a1 =  2.0e0 * y(ic12) * ratdum(ir1212)                            &
       &    - y(ihe4) * ratdum(ircag)                                     &
       &    + 0.5e0 * y(io16) * ratdum(ir1216)                            &
       &    + 3.e0 * ratdum(irg3a) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(n14)
  a1  =  y(ih1)*ratdum(ifa)*ratdum(irnpg) -y(ihe4)*ratdum(irnag)*1.5
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(o16)
  a1 =  y(ih1) * ratdum(iropg)                                      &
       &    + 1.12e0 * y(io16) * ratdum(ir1616)                           &
       &    - y(ihe4) * ratdum(iroag)                                     &
       &    + 0.5e0 * y(ic12) * ratdum(ir1216)                            &
       &    + s1 * y(io16) * 0.68e0 * ratdum(ir1616)                      &
       &    + ratdum(iroga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ne20)
  a1  = -y(ihe4) * ratdum(irneag)  + ratdum(irnega)            
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(mg24)
  a1  =  r1 * y(ihe4) * ratdum(irmgap)                              &
       &     - y(ihe4) * (ratdum(irmgag)+ratdum(irmgap))                  &
       &     + ratdum(irmgga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(si28)
  a1  =  s1 * y(ihe4) * ratdum(irsiap)                              &
       &     - y(ihe4) * (ratdum(irsiag)+ratdum(irsiap))                  &
       &     + (ratdum(irsiga)+r1 * ratdum(irsigp)) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(s32)
  a1  =  t1 * y(ihe4) * ratdum(irsap)                               &
       &     - y(ihe4) * (ratdum(irsag)+ratdum(irsap))                    &
       &     + (ratdum(irsga)+s1 * ratdum(irsgp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ar36)
  a1  =  u1 * y(ihe4) * ratdum(irarap)                              &
       &     - y(ihe4) * (ratdum(irarag)+ratdum(irarap))                  &
       &     + (ratdum(irarga)+t1 * ratdum(irargp)) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ca40)
  a1  =  v1 * y(ihe4) * ratdum(ircaap)                              &
       &     - y(ihe4) * (ratdum(ircaag)+ratdum(ircaap))                  &
       &     + (ratdum(ircaga)+u1 * ratdum(ircagp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ti44)
  a1  =  w1 * y(ihe4) * ratdum(irtiap)                              &
       &     - y(ihe4) * (ratdum(irtiag)+ratdum(irtiap))                  &
       &     + (ratdum(irtiga)+v1 * ratdum(irtigp)) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(cr48)
  a1  =  x1 * y(ihe4) * ratdum(ircrap)                              &
       &     - y(ihe4) * (ratdum(ircrag)+ratdum(ircrap))                  &
       &     + (ratdum(ircrga)+w1 * ratdum(ircrgp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(fe52)
  a1  = -y(ihe4) * y(iprot) * r7f54                                 &
       &     - y(ihe4) * (ratdum(irfeag)+r6f54)                           &
       &     + (ratdum(irfega)+x1 * ratdum(irfegp)) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(fe54)
  a1  =  yprot2 * r5f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(ni56)
  a1  =  (ratdum(irniga)+y(iprot) * r8f54)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(neut)
  a1  =  2.0e0 * y(ineut) * yprot2 * ralf2
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(he4)/d(prot)
  a1  = -y(ife52) * y(ihe4)*r7f54 + y(ini56)*r8f54                  &
       &     + 2.0e0 * yneut2 * y(iprot) * ralf2                          &
       &     + 2.0e0 * y(ife54) * y(iprot) * r5f54                        &
       &                  + y(ife54) * yprot2 * r5f54a                    &
       &                  - y(ihe4) * y(ife52) * r6f54a                   &
       &                  - y(ife52) * y(ihe4) * y(iprot) * r7f54a        &
       &                  + y(ini56) * y(iprot) * r8f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..c12 jacobian element
  !..d(c12)/d(h1)
  a1  =  -y(ic12) * ratdum(ircpg)                                   &
       &     + y(in14) * ratdum(ifa) * ratdum(irnpg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(he4)
  a1  =  -y(ic12) * ratdum(ircag)                                   &
       &      + 3.0e0 * y(ihe4) * y(ihe4) * ratdum(ir3a)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(c12)
  a1  =  -4.0e0 * y(ic12) * ratdum(ir1212)                          &
       &     - y(ihe4) * ratdum(ircag)                                    &
       &     - y(io16) * ratdum(ir1216)  - ratdum(irg3a)                  &
       &     - y(ih1) * ratdum(ircpg)          
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(n14)
  a1  =  y(ih1) * ratdum(ifa) * ratdum(irnpg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(c12)/d(o16)
  a1  = -y(ic12) * ratdum(ir1216)   + ratdum(iroga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..n14 jacobian elements
  !..d(n14)/d(h1)
  a1  =  y(ic12) * ratdum(ircpg)                                    &
       &     - y(in14) * ratdum(irnpg)                                    &
       &     + y(io16) * ratdum(iropg)         
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(n14)/d(he4)
  a1  = -y(in14) * ratdum(irnag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(n14)/d(c12)
  a1  =  y(ih1) * ratdum(ircpg)      
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(n14)/d(n14)
  a1  = -y(ih1) * ratdum(irnpg) - y(ihe4) * ratdum(irnag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(n14)/d(o16)
  a1  = y(ih1) * ratdum(iropg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..16o jacobian elements
  !..d(o16)/d(h1)
  a1  =  y(in14)*ratdum(ifg)*ratdum(irnpg) - y(io16)*ratdum(iropg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(he4)
  a1  = -y(io16) * ratdum(iroag) + y(ic12) * ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(c12)
  a1  = -y(io16) * ratdum(ir1216) + y(ihe4) * ratdum(ircag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1
  !
  !..d(o16)/d(n14)
  a1  =  y(ih1) * ratdum(ifg) * ratdum(irnpg)      
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(o16)
  a1  = -4.0e0 * y(io16) * ratdum(ir1616)                           &
       &     - y(ihe4) * ratdum(iroag)                                    &
       &     - y(ic12) * ratdum(ir1216)  - ratdum(iroga)                  &
       &     - y(ih1) * ratdum(iropg)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(o16)/d(ne20)
  a1  =  ratdum(irnega)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..20ne jacobian elements
  !..d(ne20)/d(he4)
  a1  = -y(ine20) * ratdum(irneag)                                  &
       &     + y(io16) * ratdum(iroag)                                    &
       &     + y(in14) * ratdum(irnag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(c12)
  a1  =  2.0e0 * y(ic12) * ratdum(ir1212)      
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(n14)
  a1  =  y(ihe4) * ratdum(irnag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(o16)
  a1  =  y(ihe4) * ratdum(iroag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(ne20)
  a1  = -y(ihe4) * ratdum(irneag) - ratdum(irnega)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ne20)/d(mg24)
  a1  =  ratdum(irmgga)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..24mg jacobian elements
  !..d(mg24)/d(he4)
  dum1 = ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1)
  a1  = -y(img24)*dum1 + y(ine20) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(c12)
  a1  =  0.5e0 * y(io16) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(o16)
  a1  =  0.5e0 * y(ic12) * ratdum(ir1216)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(ne20)
  a1  =  y(ihe4) * ratdum(irneag)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(mg24)
  a1  = -y(ihe4)*dum1 - ratdum(irmgga) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(mg24)/d(si28)
  a1  = ratdum(irsiga)+r1 * ratdum(irsigp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..28si jacobian elements
  !..d(si28)/d(he4)
  dum1 = ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1)
  dum2 = ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1)
  a1  =  y(img24)*dum1 - y(isi28)*dum2
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(c12)
  a1  =  0.5e0 * y(io16) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(o16)
  a1  =  1.12e0 * y(io16) * ratdum(ir1616)                          &
       &     + 0.68e0 * y(io16) * s1 * ratdum(ir1616)                     &
       &     + 0.5e0 * y(ic12) * ratdum(ir1216) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(mg24)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(si28)
  a1  = -y(ihe4)*dum2 - (r1 * ratdum(irsigp)+ratdum(irsiga))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(si28)/d(s32)
  a1  =  ratdum(irsga)+s1 * ratdum(irsgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..32s jacobian elements
  !..d(s32)/d(he4)
  dum1 = ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1)
  dum2 = ratdum(irsag)+ratdum(irsap) * (1.0e0-t1)
  a1   =  y(isi28)*dum1 - y(is32)*dum2
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(s32)/d(o16)
  a1  =  0.68e0 * y(io16)*ratdum(ir1616) * (1.0e0-s1)               &
       &     + 0.2e0 * y(io16) * ratdum(ir1616) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(s32)/d(si28)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(s32)/d(s32)
  a1  = -y(ihe4)*dum2 - (ratdum(irsga)+s1 * ratdum(irsgp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(s32)/d(ar36)
  a1  =  ratdum(irarga)+t1 * ratdum(irargp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..36ar jacobian elements
  !..d(ar36)/d(he4)
  dum1 = ratdum(irsag)+ratdum(irsap) * (1.0e0-t1)
  dum2 = ratdum(irarag)+ratdum(irarap) * (1.0e0-u1)
  a1   =  y(is32)*dum1 - y(iar36)*dum2
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ar36)/d(s32)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ar36)/d(ar36)
  a1  = -y(ihe4)*dum2 - (ratdum(irarga)+t1 * ratdum(irargp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1
  !
  !..d(ar36)/d(ca40)
  a1  =  ratdum(ircaga)+ratdum(ircagp) * u1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..40ca jacobian elements
  !..d(ca40)/d(he4)
  dum1 = ratdum(irarag)+ratdum(irarap) * (1.0e0-u1)
  dum2 = ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1)
  a1  =  y(iar36)*dum1 - y(ica40)*dum2
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ca40)/d(ar36)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ca40)/d(ca40)
  a1  = -y(ihe4)*dum2 - (ratdum(ircaga)+ratdum(ircagp) * u1)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ca40)/d(ti44)
  a1  =  ratdum(irtiga)+ratdum(irtigp) * v1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..44ti jacobian elements
  !..d(ti44)/d(he4)
  dum1 = ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1)
  dum2 = ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1)
  a1   =  y(ica40)*dum1 - y(iti44)*dum2
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ti44)/d(ca40)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ti44)/d(ti44)
  a1  = -y(ihe4)*dum2 - (ratdum(irtiga)+v1 * ratdum(irtigp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ti44)/d(cr48)
  a1  =  ratdum(ircrga)+w1 * ratdum(ircrgp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..48cr jacobian elements
  !..d(cr48)/d(he4)
  dum1 = ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1)
  dum2 = ratdum(ircrag)+ratdum(ircrap) * (1.0e0-x1)
  a1   =  y(iti44)*dum1 - y(icr48)*dum2
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(cr48)/d(ti44)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(cr48)/d(cr48)
  a1  = -y(ihe4)*dum2 - (ratdum(ircrga)+w1 * ratdum(ircrgp))
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(cr48)/d(fe52)
  a1  =  ratdum(irfega)+x1 * ratdum(irfegp)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..52fe jacobian elements
  !..d(fe52)/d(he4)
  dum1 = ratdum(ircrag)+(1.0e0-x1) * ratdum(ircrap)
  a1   =  y(icr48)*dum1 - y(ife52) * (ratdum(irfeag)+r6f54)         &
       &       - y(ife52) * y(iprot) * r7f54  
  nt   = nt + 1
  iat  = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(cr48)
  a1  =  y(ihe4)*dum1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(fe52)
  a1  = -(ratdum(irfega)+x1 * ratdum(irfegp))                       &
       &     - y(ihe4) * (ratdum(irfeag)+r6f54)                           &
       &     - y(ihe4) * y(iprot) * r7f54                                 &
       &     - yneut2 * r2f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(fe54)
  a1  =  yprot2 * r5f54 + r1f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(ni56)
  a1  =  ratdum(irniga) + y(iprot)*r8f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(neut)
  a1  = -2.0e0 * y(ife52) * y(ineut) * r2f54                        &
       &                     + y(ife54) * r1f54a                          &
       &                     - y(ife52) * yneut2 * r2f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe52)/d(prot)
  a1  =  2.0e0 * y(ife54) * y(iprot) * r5f54                        &
       &     - y(ife52)*y(ihe4)*r7f54  + y(ini56)*r8f54                   &
       &                   + y(ife54) * yprot2 * r5f54a                   &
       &                   - y(ife52) * y(ihe4) * r6f54a                  &
       &                   - y(ife52) * y(ihe4) * y(iprot) * r7f54a       &
       &                   + y(ini56) * y(iprot) * r8f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..54fe jacobian elements
  !..d(fe54)/d(he4)
  a1  =  y(ife52) * r6f54  
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe54)/d(fe52)
  a1  =  yneut2 * r2f54    + y(ihe4) * r6f54           
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe54)/d(fe54)
  a1  = -r1f54  - yprot2 * (r3f54+r5f54) 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe54)/d(ni56)
  a1  =  r4f54  + 56.0e0 * ratdum(irn56ec)/54.0e0
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1
  !
  !..d(fe54)/d(neut)
  a1  =  2.0e0 * y(ife52) * y(ineut) * r2f54                        &
       &                     - y(ife54) * r1f54a                          &
       &                     + y(ife52) * yneut2 * r2f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(fe54)/d(prot)
  a1  = -2.0e0 * y(ife54) * y(iprot) * (r3f54+r5f54)                &
       &                   - y(ife54) * yprot2 * (r3f54a+r5f54a)          &
       &                    + y(ini56) * r4f54a                           &
       &                    + y(ife52) * y(ihe4) * r6f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1



  !..56ni jacobian elements
  !..d(ni56)/d(he4)
  a1  =  y(ife52) * (ratdum(irfeag)+y(iprot) * r7f54)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(fe52)
  a1  =  y(ihe4) * (ratdum(irfeag)+y(iprot) * r7f54)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(fe54)
  a1  =  yprot2 * r3f54  
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(ni56)
  a1  = -(ratdum(irniga) + r4f54 + y(iprot)*r8f54) - ratdum(irn56ec)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(ni56)/d(prot)
  a1  =  y(ife52) * y(ihe4) * r7f54  - y(ini56)*r8f54               &
       &     + 2.0e0 * y(ife54) * y(iprot) * r3f54                        &
       &                   + y(ife54) * yprot2 * r3f54a                   &
       &                   - y(ini56) * (r4f54a + y(iprot)*r8f54a)        &
       &                   + y(ife52) * y(ihe4) * y(iprot) * r7f54a
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..neutron jacobian elements
  !..d(neut)/d(he4)
  a1  =  2.0e0 * ralf1
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(fe52)
  a1  = -2.0e0 * yneut2 * r2f54 
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(fe54)
  a1  =  2.0e0 * r1f54        
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(neut)
  a1  =  -4.0e0 * y(ife52) * y(ineut) * r2f54                       &
       &       - 4.0e0 * yprot2 * y(ineut) * ralf2                        &
       &       + 2.0e0 * y(ife54) * r1f54a                                &
       &       - 2.0e0 * y(ife52) * yneut2 * r2f54a                       &
       &       - ratdum(irnep)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(prot)
  a1  = -4.0e0 * y(iprot) * yneut2 * ralf2                          &
       &                   + ratdum(irpen)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..proton jacobian elements
  !..d(prot)/d(he4)
  a1  =  2.0e0 * ralf1   + 2.0e0 * y(ife52) * r6f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(fe52)
  a1  =  2.0e0 * y(ihe4) * r6f54
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(fe545)
  a1  = -2.0e0 * yprot2 * (r3f54+r5f54)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(ni56)
  a1  =  2.0e0 * r4f54           
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(neut)
  a1  = -4.0e0 * yprot2 * y(ineut) * ralf2                          &
       &                    + ratdum(irnep)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1

  !..d(neut)/d(prot)
  a1  = -4.0e0 * y(iprot) * yneut2 * ralf2                          &
       &     - 4.0e0 * y(ife54) * y(iprot) * (r3f54+r5f54)                &
       &     - 2.0e0 * y(ife54) * yprot2 * (r3f54a+r5f54a)                &
       &     + 2.0e0 * y(ini56) * r4f54a                                  &
       &     + 2.0e0 * y(ife52) * y(ihe4) * r6f54a                        &
       &     - ratdum(irpen)
  nt  = nt + 1
  iat = eloc(nt)
  dfdy(iat) = dfdy(iat) + a1


  !..bullet check the counting
  if (nt .ne. nterms) then
     write(6,*) 'nt =',nt,'  nterms =',nterms
     write(6,*) 'error in routine bn_networkSparseJakob: nt .ne. nterms'
     call Driver_abortFlash('ERROR in bn_networkSparseJakob: nt /= nterms')
  end if
  return
end subroutine bn_networkSparseJakob






