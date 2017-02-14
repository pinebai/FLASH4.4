!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_network
!!
!! NAME
!! 
!! bn_network
!!
!! SYNOPSIS
!!
!! bn_network(real, intent(IN)                :: tt,
!!            real, intent(IN), dimension(:)  :: y,
!!            real, intent(OUT), dimension(:) :: dydt
!!
!! DESCRIPTION
!!
!!  this routine sets up the system of ode's for the aprox13 nuclear reactions.
!!  this is an alpha chain + heavy ion network with (a,p)(p,g) links
!!  
!!  isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
!!            ar36, ca40, ti44, cr48, fe52, ni56
!!
!! ARGUMENTS
!!
!!   tt  -- never used, so don't know what it does
!!   y --
!!   dydt --
!!
!!***

subroutine bn_network(tt,y,dydt)   

!! Sorry, there are just toooo freakin many of them to include by name
  use Burn_data  
  use bn_dataAprox13

  implicit none

#include "constants.h"
#include "Flash.h"


  ! branching ratios
 
!!  declare arguments and local variables
  real, intent(IN)  :: tt,y(*)
  real, intent(OUT) :: dydt(*)

  real             :: r1,s1,t1,u1,v1,w1,x1,y1


  !!  branching ratios for various dummy proton links
  r1     = ratdum(iralpa)/(ratdum(iralpa) + ratdum(iralpg))
  s1     = ratdum(irppa)/(ratdum(irppa) + ratdum(irppg))
  t1     = ratdum(irclpa)/(ratdum(irclpa) + ratdum(irclpg))
  u1     = ratdum(irkpa)/(ratdum(irkpa) + ratdum(irkpg))
  v1     = ratdum(irscpa)/(ratdum(irscpa) + ratdum(irscpg))
  w1     = ratdum(irvpa)/(ratdum(irvpa) + ratdum(irvpg))
  x1     = ratdum(irmnpa)/(ratdum(irmnpa) + ratdum(irmnpg))
  y1     = ratdum(ircopa)/(ratdum(ircopa) + ratdum(ircopg))


  !!  set up the system of odes: 
  !!  he4 reactions
  !!  heavy ion reactions
  dydt(ihe4) =  y(ic12) * y(ic12) * ratdum(ir1212)                  &
       &            + 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)          &
       &            + 0.56e0 * y(io16) * y(io16) * ratdum(ir1616)         &
       &            + 0.34e0 * s1 * y(io16) * y(io16) * ratdum(ir1616)

  !!  (a,g) and (g,a) reactions
  dydt(ihe4) =  dydt(ihe4)                                          &
       &            - 3.0e0 * y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)  &
       &            + 3.0e0 * y(ic12) * ratdum(irg3a)                     &
       &            - y(ihe4)  * y(ic12) * ratdum(ircag)                  &
       &            + y(io16)  * ratdum(iroga)                            &
       &            - y(ihe4)  * y(io16) * ratdum(iroag)                  &
       &            + y(ine20) * ratdum(irnega)                           &
       &            - y(ihe4)  * y(ine20) * ratdum(irneag)                &
       &            + y(img24) * ratdum(irmgga)                           &
       &            - y(ihe4)  * y(img24)* ratdum(irmgag)                 &
       &            + y(isi28) * ratdum(irsiga)                           &
       &            - y(ihe4)  * y(isi28)*ratdum(irsiag)                  &
       &            + y(is32)  * ratdum(irsga)

  dydt(ihe4) =  dydt(ihe4)                                          &
       &            - y(ihe4)  * y(is32) * ratdum(irsag)                  &
       &            + y(iar36) * ratdum(irarga)                           &
       &            - y(ihe4)  * y(iar36)*ratdum(irarag)                  &
       &            + y(ica40) * ratdum(ircaga)                           &
       &            - y(ihe4)  * y(ica40)*ratdum(ircaag)                  &
       &            + y(iti44) * ratdum(irtiga)                           &
       &            - y(ihe4)  * y(iti44)*ratdum(irtiag)                  &
       &            + y(icr48) * ratdum(ircrga)                           &
       &            - y(ihe4)  * y(icr48)*ratdum(ircrag)                  &
       &            + y(ife52) * ratdum(irfega)                           &
       &            - y(ihe4)  * y(ife52) * ratdum(irfeag)                &
       &            + y(ini56) * ratdum(irniga)


  !!  (a,p)(p,g) and (g,p)(p,a) reactions
  dydt(ihe4) =  dydt(ihe4)                                          &
       &            - y(ihe4)  * y(img24) * ratdum(irmgap) * (1.0e0-r1)   &
       &            + y(isi28) * ratdum(irsigp) * r1                      &
       &            - y(ihe4)  * y(isi28) * ratdum(irsiap) * (1.0e0-s1)   &
       &            + y(is32)  * ratdum(irsgp) * s1                       &
       &            - y(ihe4)  * y(is32) * ratdum(irsap) * (1.0e0-t1)     &
       &            + y(iar36) * ratdum(irargp) * t1                      &
       &            - y(ihe4)  * y(iar36) * ratdum(irarap) * (1.0e0-u1)   &
       &            + y(ica40) * ratdum(ircagp) * u1                      &
       &            - y(ihe4)  * y(ica40) * ratdum(ircaap) * (1.0e0-v1)   &
       &            + y(iti44) * ratdum(irtigp) * v1                      &
       &            - y(ihe4)  * y(iti44) * ratdum(irtiap) * (1.0e0-w1)   &
       &            + y(icr48) * ratdum(ircrgp) * w1                      &
       &            - y(ihe4)  * y(icr48) * ratdum(ircrap) * (1.0e0-x1)   &
       &            + y(ife52) * ratdum(irfegp) * x1                      &
       &            - y(ihe4)  * y(ife52) * ratdum(irfeap) * (1.0e0-y1)   &
       &            + y(ini56) * ratdum(irnigp) * y1


  !!  c12 reactions   
  dydt(ic12) = -2.0e0 * y(ic12) * y(ic12) * ratdum(ir1212)          &
       &            - y(ic12) * y(io16) * ratdum(ir1216)                  &
       &            + y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)          &
       &            - y(ic12) * ratdum(irg3a)                             &
       &            - y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            + y(io16) * ratdum(iroga)

  !!  16o reactions 
  dydt(io16) = -y(ic12) * y(io16) * ratdum(ir1216)                  &
       &            - 2.0e0 * y(io16) * y(io16) * ratdum(ir1616)          &
       &            + y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            - y(io16) * y(ihe4) * ratdum(iroag)                   &
       &            - y(io16) * ratdum(iroga)                             &
       &            + y(ine20) * ratdum(irnega)

  !!  20ne reactions 
  dydt(ine20) =  y(ic12) * y(ic12) * ratdum(ir1212)                 &
       &             + y(io16) * y(ihe4) * ratdum(iroag)                  &
       &             - y(ine20) * y(ihe4) * ratdum(irneag)                &
       &             - y(ine20) * ratdum(irnega)                          &
       &             + y(img24) * ratdum(irmgga)            

  !!  24mg reactions
  dydt(img24)  = 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)         &
       &              + y(ine20) * y(ihe4) * ratdum(irneag)               &
       &              - y(img24) * y(ihe4) * ratdum(irmgag)               &
       &              - y(img24) * ratdum(irmgga)                         &
       &              + y(isi28) * ratdum(irsiga)                         &
       &              - y(img24) *  y(ihe4) *  ratdum(irmgap) * (1.0e0-r1)&
       &              + y(isi28) * r1 * ratdum(irsigp)

  !!  28si reactions
  dydt(isi28)  =  0.5e0 * y(ic12) * y(io16) * ratdum(ir1216) &
       &              + 0.56e0 * y(io16) * y(io16) * ratdum(ir1616)&
       &              + 0.34e0 * y(io16) * y(io16) * s1 * ratdum(ir1616)&
       &              + y(img24) * y(ihe4) * ratdum(irmgag) &
       &              - y(isi28) * y(ihe4) * ratdum(irsiag) &
       &              - y(isi28) * ratdum(irsiga)&
       &              + y(is32)  * ratdum(irsga) &
       &              + y(img24) * y(ihe4) * ratdum(irmgap) * (1.0e0-r1)&
       &              - y(isi28) * r1 * ratdum(irsigp) &
       &              - y(isi28) * y(ihe4) * ratdum(irsiap) * (1.0e0-s1)&
       &              + y(is32)  * s1 * ratdum(irsgp)

  !!  32s reactions
  dydt(is32) = + 0.34e0 * y(io16)*y(io16)*ratdum(ir1616)*(1.0e0-s1) &
       &             + 0.1e0 * y(io16) * y(io16) * ratdum(ir1616)  &
       &             +  y(isi28) * y(ihe4) * ratdum(irsiag)  &
       &             - y(is32) * y(ihe4) * ratdum(irsag)  &
       &             - y(is32) * ratdum(irsga)  &
       &             + y(iar36) * ratdum(irarga) &
       &             + y(isi28) * y(ihe4) * ratdum(irsiap) * (1.0e0-s1) &
       &             - y(is32) * s1 * ratdum(irsgp) &
       &             - y(is32) * y(ihe4) * ratdum(irsap) * (1.0e0-t1) &
       &             + y(iar36) * t1 * ratdum(irargp)


  !!  36ar reactions
  dydt(iar36) =  y(is32)  * y(ihe4) * ratdum(irsag) &
       &             - y(iar36) * y(ihe4) * ratdum(irarag) &
       &             - y(iar36) * ratdum(irarga)  &
       &             + y(ica40) * ratdum(ircaga) &
       &             + y(is32)  * y(ihe4) * ratdum(irsap) * (1.0e0-t1)  &
       &             - y(iar36) * t1 * ratdum(irargp) &
       &             - y(iar36) * y(ihe4) * ratdum(irarap) * (1.0e0-u1) &
       &             + y(ica40) * ratdum(ircagp) * u1


  !!  40ca reactions
  dydt(ica40) =  y(iar36) * y(ihe4) * ratdum(irarag) &
       &             - y(ica40) * y(ihe4) * ratdum(ircaag) &
       &             - y(ica40) * ratdum(ircaga)  &
       &             + y(iti44) * ratdum(irtiga)  &
       &             + y(iar36) * y(ihe4) * ratdum(irarap) * (1.0e0-u1) &
       &             - y(ica40) * ratdum(ircagp) * u1 &
       &             - y(ica40) * y(ihe4) * ratdum(ircaap) * (1.0e0-v1) &
       &             + y(iti44) * ratdum(irtigp) * v1

  !!  44ti reactions
  dydt(iti44) =  y(ica40) * y(ihe4) * ratdum(ircaag) &
       &             - y(iti44) * y(ihe4) * ratdum(irtiag) &
       &             - y(iti44) * ratdum(irtiga)  &
       &             + y(icr48) * ratdum(ircrga) &
       &             + y(ica40) * y(ihe4) * ratdum(ircaap)*(1.0e0-v1) &
       &             - y(iti44) * v1 * ratdum(irtigp) &
       &             - y(iti44) * y(ihe4) * ratdum(irtiap) * (1.0e0-w1) &
       &             + y(icr48) * w1 * ratdum(ircrgp)

  !!  48cr reactions
  dydt(icr48) =  y(iti44) * y(ihe4) * ratdum(irtiag) &
       &             - y(icr48) * y(ihe4) * ratdum(ircrag) &
       &             - y(icr48) * ratdum(ircrga)  &
       &             + y(ife52) * ratdum(irfega)  &
       &             + y(iti44) * y(ihe4) * ratdum(irtiap)*(1.0e0-w1) &
       &             - y(icr48) * w1 * ratdum(ircrgp) &
       &             - y(icr48) * y(ihe4) * ratdum(ircrap) * (1.0e0-x1) &
       &             + y(ife52) * x1 * ratdum(irfegp)

  !!  52fe reactions
  dydt(ife52) =  y(icr48) * y(ihe4) * ratdum(ircrag) &
       &             - y(ife52) * y(ihe4) * ratdum(irfeag) &
       &             - y(ife52) * ratdum(irfega) &
       &             + y(ini56) * ratdum(irniga) &
       &             + y(icr48) * y(ihe4) * ratdum(ircrap) * (1.0e0-x1)  &
       &             - y(ife52) * x1 * ratdum(irfegp) &
       &             - y(ife52) * y(ihe4) * ratdum(irfeap) * (1.0e0-y1) &
       &             + y(ini56) * y1 * ratdum(irnigp)

  !! 56ni reactions
  dydt(ini56) =  y(ife52) * y(ihe4) * ratdum(irfeag) &
       &             - y(ini56) * ratdum(irniga) &
       &             + y(ife52) * y(ihe4) * ratdum(irfeap) * (1.0e0-y1)   &
       &             - y(ini56) * y1 * ratdum(irnigp) 

  return

end   subroutine bn_network

