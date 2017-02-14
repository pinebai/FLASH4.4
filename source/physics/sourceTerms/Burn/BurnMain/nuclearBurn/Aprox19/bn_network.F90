!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_network
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
!!  Aprox19 was the third network installed in flash 1.5
!!  it is like the original iso13 network, but adds protons and he3
!!  for pp burning, nitrogen to do cno and hot cno (beta limited)
!!  burning, fe54 and photodisintegration protons and neutrons
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

  use Burn_dataEOS, ONLY:  btemp
  use Burn_data
  use bn_dataAprox19
 
  implicit none

  !..this routine sets up the system of odes for the network nuclear reactions.
  !..this is an alpha chain with (a,p)(p,g) links, heavy ion, pp, cno,
  !..and photodisintegration network.
  !..
  !..isotopes: h1, he3, he4, c12, n14, o16, ne20, mg24, si28, s32, ar36,
  !..          ca40, ti44, cr48, fe52, fe54, ni56, neut, prot


!!  declare arguments and local variables
  real, intent(IN)  :: tt,y(*)
  real, intent(OUT) :: dydt(*)

  real         yneut2,yprot2,xx,den1,den2,                &
       &       r1,s1,t1,u1,v1,w1,x1,ralf1,ralf2,                          &
       &       r1f54,r2f54,r3f54,r4f54,r5f54,r6f54,r7f54,r8f54,           &
       &       yy,den1a,den2a,r1f54a,r2f54a,                              &
       &       r3f54a,r4f54a,r5f54a,r6f54a,r7f54a,r8f54a



  !..branching ratios for dummy proton links
  !..and special combined rates for alpha photodisintegration and fe54
  real :: entropy, dst, dsd
  yneut2 = y(ineut) * y(ineut)
  yprot2 = y(iprot) * y(iprot)
  r1     = ratdum(iralpa)/(ratdum(iralpa) + ratdum(iralpg))
  s1     = ratdum(irppa)/(ratdum(irppa) + ratdum(irppg))
  t1     = ratdum(irclpa)/(ratdum(irclpa) + ratdum(irclpg))
  u1     = ratdum(irkpa)/(ratdum(irkpa) + ratdum(irkpg))
  v1     = ratdum(irscpa)/(ratdum(irscpa) + ratdum(irscpg))
  w1     = ratdum(irvpa)/(ratdum(irvpa) + ratdum(irvpg))
  x1     = ratdum(irmnpa)/(ratdum(irmnpa) + ratdum(irmnpg))
  xx     = ratdum(irhegp)*ratdum(irdgn) +                           &
       &         y(ineut)*ratdum(irheng)*ratdum(irdgn) +                  &
       &         y(ineut)*y(iprot)*ratdum(irheng)*ratdum(irdpg)
  ralf1  = 0.0e0
  ralf2  = 0.0e0
  if (xx .ge. 1.0e-150) then
     ralf1 = ratdum(irhegn)*ratdum(irhegp)*ratdum(irdgn)/xx
     ralf2 = ratdum(irheng)*ratdum(irdpg)*ratdum(irhng)/xx
  end if


  !..don't consider 54fe photodisintegration if temperature is too low or there 
  !..is still free hydrogen around.
  xx  = ratdum(ircogp) + y(iprot)*(ratdum(ircopg) + ratdum(ircopa))
  if (xx .gt. 1.0e-99) then
     den1  = 1.0/xx
  else
     den1  = 0.0e0
  end if
  if (btemp/1.0e9 .lt. 1.5 .or. y(ih1) .gt. 1.0e-5) then
     den1  = 0.0e0
  end if

  yy  = ratdum(ir53gn) + y(ineut)*ratdum(ir53ng)
  if (yy .gt. 1.0e-99) then
     den2  = 1.0/yy
  else
     den2  = 0.0e0
  end if
  if (btemp/1.0e9 .lt. 1.5 .or. y(ih1) .gt. 1.0e-5) then
     den2  = 0.0e0
  end if

  !       den1 = 0.0e0
  !       den2 = 0.0e0


  r1f54  = ratdum(ir54gn)*ratdum(ir53gn)*den2
  r2f54  = ratdum(ir52ng)*ratdum(ir53ng)*den2
  r3f54  = ratdum(irfepg)*ratdum(ircopg)*den1
  r4f54  = ratdum(irnigp)*ratdum(ircogp)*den1
  r5f54  = ratdum(irfepg)*ratdum(ircopa)*den1
  r6f54  = ratdum(irfeap)*ratdum(ircogp)*den1
  r7f54  = ratdum(irfeap)*ratdum(ircopg)*den1
  r8f54  = ratdum(irnigp)*ratdum(ircopa)*den1



  !..beta limit the he3+he4, 14n(pg), and 16o(pg) rates. assume steady state 
  !..abundances of 14o and 15o in beta limited cno cycle.

  if (y(ihe4) .gt. 0.0)                                             &
       &   ratdum(ir34)  = min(ratdum(ir34),0.896e0/y(ihe4))
  if (y(ih1)  .gt. 0.0)                                             &
       &   ratdum(irnpg) = min(ratdum(irnpg),5.68e-03/(y(ih1)*1.57e0))
  if (y(ih1)  .gt. 0.0)                                             &
       &   ratdum(iropg) = min(ratdum(iropg),0.0105e0/y(ih1))




  !..set up the system of ode's : 
  !..hydrogen reactions   
  dydt(ih1) =  -3.0e0 * y(ih1) * y(ih1) * ratdum(irpp)              &
       &            + 2.0e0 * y(ihe3) * y(ihe3) * ratdum(ir33)            &
       &            - y(ihe3) * y(ihe4) * ratdum(ir34)                    &
       &            - 2.0e0 * y(ic12) * y(ih1) * ratdum(ircpg)            &
       &            - 2.0e0 * y(in14) * y(ih1) * ratdum(irnpg)            &
       &            - 2.0e0 * y(io16) * y(ih1) * ratdum(iropg)            &
       &            - 3.0e0 * y(ih1) * ratdum(irpen)


  !..he3 reactions 
  dydt(ihe3) =  y(ih1) * y(ih1) * ratdum(irpp)                      &
       &            - 2.0e0 * y(ihe3) * y(ihe3) * ratdum(ir33)            &
       &            - y(ihe3) * y(ihe4) * ratdum(ir34)                    &
       &            + y(ih1) * ratdum(irpen)


  !..he4 reactions
  dydt(ihe4) =  y(ihe3) * y(ihe3) * ratdum(ir33)                    &
       &            + y(ihe3) * y(ihe4) * ratdum(ir34)                    &
       &            + y(in14) * y(ih1) * ratdum(ifa) * ratdum(irnpg)      &
       &            + y(io16) * y(ih1) * ratdum(iropg)                    &
       &            - y(ihe4) * y(in14) * ratdum(irnag) * 1.5             &
       &            + y(ic12) * y(ic12) * ratdum(ir1212)                  &
       &            + 0.56 * y(io16) * y(io16) * ratdum(ir1616)           &
       &            + r1 * y(ihe4) * y(img24) * ratdum(irmgap)            &
       &            + s1 * y(ihe4) * y(isi28) * ratdum(irsiap)            &
       &            + t1 * y(ihe4) * y(is32) * ratdum(irsap)              &
       &            + u1 * y(ihe4) * y(iar36) * ratdum(irarap)            &
       &            - y(ihe4) * y(io16) * ratdum(iroag)                   &
       &            - y(ihe4) * y(ine20) * ratdum(irneag)                 &
       &            - y(ihe4) * y(ic12) * ratdum(ircag)                   &
       &            - y(ihe4) * y(img24) * (ratdum(irmgag)+ratdum(irmgap))&
       &            - y(ihe4) * y(isi28) * (ratdum(irsiag)+ratdum(irsiap))
  dydt(ihe4) =  dydt(ihe4)                                          &
       &            - y(ihe4) * y(is32) * (ratdum(irsag)+ratdum(irsap))   &
       &            - y(ihe4) * y(iar36) * (ratdum(irarag)+ratdum(irarap))&
       &            + 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)          &
       &            + v1 * y(ihe4) * y(ica40) * ratdum(ircaap)            &
       &            + w1 * y(ihe4) * y(iti44) * ratdum(irtiap)            &
       &            + x1 * y(ihe4) * y(icr48) * ratdum(ircrap)            &
       &            - y(ife52) * y(ihe4) * y(iprot) * r7f54               &
       &            - y(ihe4) * y(ica40) * (ratdum(ircaag)+ratdum(ircaap))&
       &            - y(ihe4) * y(iti44) * (ratdum(irtiag)+ratdum(irtiap))&
       &            - y(ihe4) * y(icr48) * (ratdum(ircrag)+ratdum(ircrap))&
       &            - y(ihe4) * y(ife52) * (ratdum(irfeag)+r6f54)         &
       &            - 3.0e0 * y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)  &
       &            + s1 * y(io16) * y(io16) * 0.34 * ratdum(ir1616)      &
       &            + 3.e0 * y(ic12) * ratdum(irg3a)                      &
       &            + y(io16) * ratdum(iroga)
  dydt(ihe4) =  dydt(ihe4)                                          &
       &            + y(ine20) * ratdum(irnega)                           &
       &            + y(img24) * ratdum(irmgga)                           &
       &            + y(isi28) * (ratdum(irsiga) + r1 * ratdum(irsigp))   &
       &            + y(is32) * (ratdum(irsga) + s1 * ratdum(irsgp))      &
       &            + y(iar36) * (ratdum(irarga) + t1 * ratdum(irargp))   &
       &            + y(ica40) * (ratdum(ircaga) + u1 * ratdum(ircagp))   &
       &            + y(iti44) * (ratdum(irtiga) + v1 * ratdum(irtigp))   &
       &            + y(icr48) * (ratdum(ircrga) + w1 * ratdum(ircrgp))   &
       &            + y(ife52) * (ratdum(irfega) + x1 * ratdum(irfegp))   &
       &            + y(ini56) * (ratdum(irniga) + y(iprot) * r8f54)      &
       &            - y(ihe4) * ralf1                                     &
       &            + yneut2 * yprot2 * ralf2                             &
       &            + y(ife54) * yprot2 * r5f54


  !..c12 reactions   
  dydt(ic12) = -2.0e0 * y(ic12) * y(ic12) * ratdum(ir1212)          &
       &            - y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            - y(ic12) * y(io16) * ratdum(ir1216)                  &
       &            + y(ihe4) * y(ihe4) * y(ihe4) * ratdum(ir3a)          &
       &            - y(ic12) * ratdum(irg3a)                             &
       &            + y(io16) * ratdum(iroga)                             &
       &            - y(ic12) * y(ih1) * ratdum(ircpg)                    &
       &            + y(in14) * y(ih1) * ratdum(ifa) * ratdum(irnpg)


  !..n14 reactions
  dydt(in14) =  y(ic12) * y(ih1) * ratdum(ircpg)                    &
       &            - y(in14) * y(ih1) * ratdum(irnpg)                    &
       &            + y(io16) * y(ih1) * ratdum(iropg)                    &
       &            - y(ihe4) * y(in14) * ratdum(irnag)


  !..16o reactions 
  dydt(io16) = -2.0e0 * y(io16) * y(io16) * ratdum(ir1616)          &
       &            - y(io16) * y(ihe4) * ratdum(iroag)                   &
       &            - y(ic12) * y(io16) * ratdum(ir1216)                  &
       &            + y(ic12) * y(ihe4) * ratdum(ircag)                   &
       &            - y(io16) * ratdum(iroga)                             &
       &            + y(ine20) * ratdum(irnega)                           &
       &            + y(in14) * y(ih1) * ratdum(ifg) * ratdum(irnpg)      &
       &            - y(io16) * y(ih1) * ratdum(iropg)


  !..20ne reactions 
  dydt(ine20) =  y(ic12) * y(ic12) * ratdum(ir1212)                 &
       &             - y(ine20) * y(ihe4) * ratdum(irneag)                &
       &             + y(io16) * y(ihe4) * ratdum(iroag)                  &
       &             - y(ine20) * ratdum(irnega)                          &
       &             + y(img24) * ratdum(irmgga)                          &
       &             + y(in14) * y(ihe4) * ratdum(irnag)



  !..24mg reactions
  dydt(img24)  = -y(img24)  *  y(ihe4)  *                           &
       &                (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1))      &
       &              + 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)        &
       &              + y(ine20) * y(ihe4) * ratdum(irneag)               &
       &              - y(img24) * ratdum(irmgga)                         &
       &              + y(isi28) * (ratdum(irsiga)                        &
       &              + r1 * ratdum(irsigp))


  !..28si reactions
  dydt(isi28)  =  y(img24)  *  y(ihe4)  *                           &
       &                (ratdum(irmgag)+ratdum(irmgap) * (1.0e0-r1))      &
       &              + 0.56e0 * y(io16) * y(io16) * ratdum(ir1616)       &
       &              - y(isi28) * y(ihe4) * (ratdum(irsiag)              &
       &              + ratdum(irsiap) * (1.0e0-s1))                      &
       &              + 0.34e0 * y(io16) * y(io16) * s1 * ratdum(ir1616)  &
       &              + 0.5e0 * y(ic12) * y(io16) * ratdum(ir1216)        &
       &              - y(isi28) * (r1 * ratdum(irsigp)+ratdum(irsiga))   &
       &              + y(is32) * (ratdum(irsga)+s1 * ratdum(irsgp))



  !..32s reactions
  dydt(is32) =  y(isi28)  *  y(ihe4)  *                             &
       &              (ratdum(irsiag)+ratdum(irsiap) * (1.0e0-s1))        &
       &             + 0.34 * y(io16)*y(io16)*ratdum(ir1616) * (1.0e0-s1) &
       &             - y(is32) * y(ihe4) * (ratdum(irsag)                 &
       &             + ratdum(irsap) * (1.0e0-t1))                        &
       &             + 0.1e0 * y(io16) * y(io16) * ratdum(ir1616)         &
       &             - y(is32) * (ratdum(irsga)+s1 * ratdum(irsgp))       &
       &             + y(iar36) * (ratdum(irarga)+t1 * ratdum(irargp))



  !..36ar reactions
  dydt(iar36) =  y(is32)  *  y(ihe4)  *                             &
       &               (ratdum(irsag)+ratdum(irsap) * (1.0e0-t1))         &
       &             - y(iar36)  *  y(ihe4)  *                            &
       &               (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1))       &
       &             - y(iar36) * (ratdum(irarga)+t1 * ratdum(irargp))    &
       &             + y(ica40) * (ratdum(ircaga)+ratdum(ircagp) * u1)



  !..40ca reactions
  dydt(ica40) =  y(iar36)  *  y(ihe4)  *                            &
       &               (ratdum(irarag)+ratdum(irarap) * (1.0e0-u1))       &
       &             - y(ica40)  *  y(ihe4)  *                            &
       &               (ratdum(ircaag)+ratdum(ircaap) * (1.0e0-v1))       &
       &             - y(ica40) * (ratdum(ircaga)+ratdum(ircagp) * u1)    &
       &             + y(iti44) * (ratdum(irtiga)+ratdum(irtigp) * v1)



  !..44ti reactions
  dydt(iti44) =  y(ica40) * y(ihe4) *                               &
       &               (ratdum(ircaag)+ratdum(ircaap)*(1.0e0-v1))         &
       &             - y(iti44) * y(ihe4) *                               &
       &               (ratdum(irtiag)+ratdum(irtiap) * (1.0e0-w1))       &
       &             - y(iti44) * (ratdum(irtiga)+v1 * ratdum(irtigp))    &
       &             + y(icr48) * (ratdum(ircrga)+w1 * ratdum(ircrgp))



  !..48cr reactions
  dydt(icr48) =  y(iti44) * y(ihe4) *                               &
       &               (ratdum(irtiag)+ratdum(irtiap)*(1.0e0-w1))         &
       &             - y(icr48) * y(ihe4) *                               &
       &               (ratdum(ircrag)+ratdum(ircrap) * (1.0e0-x1))       &
       &             - y(icr48) * (ratdum(ircrga)+w1 * ratdum(ircrgp))    &
       &             + y(ife52) * (ratdum(irfega)+x1 * ratdum(irfegp))



  !..52fe reactions
  dydt(ife52) =  y(icr48) * y(ihe4) *                               &
       &               (ratdum(ircrag)+(1.0e0-x1) * ratdum(ircrap))       &
       &             - y(ife52) * (ratdum(irfega)+x1 * ratdum(irfegp))    &
       &             - y(ife52) * y(ihe4) * (ratdum(irfeag)+r6f54)        &
       &             + y(ini56) * ratdum(irniga)                          &
       &             + y(ife54) * yprot2 * r5f54                          &
       &             - y(ife52) * y(ihe4) * y(iprot) * r7f54              &
       &             + y(ini56) * y(iprot) * r8f54                        &
       &             - y(ife52) * yneut2 * r2f54                          &
       &             + y(ife54) * r1f54



  !..54fe reactions; a double neutron/proton link to fe52/ni56
  dydt(ife54) = -y(ife54) * r1f54                                   &
       &             + y(ife52) * yneut2 * r2f54                          &
       &             - y(ife54) * yprot2 * (r3f54+r5f54)                  &
       &             + y(ini56) * r4f54                                   &
       &             + y(ife52) * y(ihe4) * r6f54                         &
       &             + 56.0e0 * ratdum(irn56ec) * y(ini56)/54.0e0



  !..56ni reactions
  dydt(ini56) =  y(ife52) * y(ihe4)*(ratdum(irfeag)+y(iprot)*r7f54) &
       &            - y(ini56) * (ratdum(irniga) + r4f54 + y(iprot)*r8f54)&
       &             + y(ife54) * yprot2 * r3f54                          &
       &             - ratdum(irn56ec) * y(ini56)



  !..neutrons from alpha photodisintegration,fe54
  dydt(ineut)  = -2.0e0 * y(ife52) * yneut2 * r2f54                 &
       &              - 2.0e0 * yprot2 * yneut2 * ralf2                   &
       &              - y(ineut) * ratdum(irnep)                          &
       &              + y(iprot) * ratdum(irpen)                          &
       &              + 2.0e0 * y(ife54) * r1f54                          &
       &              + 2.0e0 * y(ihe4) * ralf1


  !..protons from alpha photodisintegration,fe54,and weak interactions
  dydt(iprot) = -2.0e0 * yprot2 * yneut2 * ralf2                    &
       &             + 2.0e0 * y(ihe4) * ralf1                            &
       &             - y(iprot) * ratdum(irpen)                           &
       &             + y(ineut) * ratdum(irnep)                           &
       &             - 2.0e0 * y(ife54) * yprot2 * (r3f54+r5f54)          &
       &             + 2.0e0 * y(ini56) * r4f54                           &
       &             + 2.0e0 * y(ife52) * y(ihe4) * r6f54

  return
end subroutine bn_network




