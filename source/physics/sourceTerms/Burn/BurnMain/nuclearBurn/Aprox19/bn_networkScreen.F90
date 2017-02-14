!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkScreen
!!
!! NAME
!! 
!! bn_networkScreen
!!
!! SYNOPSIS
!!
!! call bn_networkScreen(real intent(IN) :: y)
!!
!! DESCRIPTION
!!
!! routine bn_networkScreen applies screening corrections to the raw reaction rates
!!
!!  this routine computes the screening factors
!!  and applies them to the raw reaction rates,
!!  producing the final reaction rates used by the
!!  right hand sides and jacobian matrix elements
!!  input y is ymass
!!
!! ARGUMENTS
!!
!!  y :  ymass
!!  
!! NOTES
!!
!! Most nuclear networks do not use the routine bn_networkWeak, called directly
!!   before bn_networkScreen.  So Aprox19 network simply calls it immediately from here.
!!
!!***
subroutine bn_networkScreen(y)
   
#include "Flash.h"

  use bn_interface, ONLY:   bn_screen4

  ! nothing seems to be passed with eos_common.fh
  !      use Burn_dataEOS
  use Burn_data
  use bn_dataAprox19

  implicit none

  !..this routine computes the screening factors
  !..and applies them to the raw reaction rates,
  !..producing the final reaction rates used by the
  !..right hand sides and jacobian matrix elements

  !!  declare
  real, intent(IN) :: y(NSPECIES)

  !! local declarations
  integer          i,j,k,jscr,screen_init,screen_on
  parameter        (screen_on = 1)
  real             sc1a,sc2a,sc3a,                   &
       &                 abar,zbar,z2bar,ytot1,zbarxx,z2barxx
  real :: entropy, dst, dsd

  !..initialize
  data             screen_init/1/

!! Call routine that is usually missing from bn_burner call
  call bn_networkWeak(y)


  !..if screening is off
  if (screen_on .eq. 0) then
     do i=1,nrat
        ratdum(i) = ratraw(i)
        scfac(i)  = 1.0e0
     end do
     return
  end if


  !..screening is on
  !..with the passed composition, compute abar,zbar and other variables
  zbarxx  = 0.0e0
  z2barxx = 0.0e0
  ytot1   = 0.0e0
  do i=1,NSPECIES
     ytot1    = ytot1 + y(i)
     zbarxx   = zbarxx + zion(i) * y(i)
     z2barxx  = z2barxx + zion(i) * zion(i) * y(i)
  enddo
  abar   = 1.0e0/ytot1
  zbar   = zbarxx * abar
  z2bar  = z2barxx * abar


  !..pp
  jscr = 1
  call bn_screen4(zbar,abar,z2bar,                          &
       &             zion(ih1),aion(ih1),zion(ih1),aion(ih1),  &
       &             jscr,screen_init,sc1a)

  ratdum(irpp)   = ratraw(irpp) * sc1a
  scfac(irpp)    = sc1a

  !..he3 + he3
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ihe3),aion(ihe3),zion(ihe3),aion(ihe3), &
       &             jscr,screen_init,sc1a)

  ratdum(ir33)   = ratraw(ir33) * sc1a
  scfac(ir33)    = sc1a

  !..he3 + he4
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                              &
       &             zion(ihe3),aion(ihe3),zion(ihe4),aion(ihe4),  &
       &             jscr,screen_init,sc1a)

  ratdum(ir34)   = ratraw(ir34) * sc1a
  scfac(ir34)    = sc1a


  !..c12 + p 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                              &
       &             zion(ic12),aion(ic12),zion(ih1),aion(ih1),    &
       &             jscr,screen_init,sc1a)

  ratdum(ircpg)  = ratraw(ircpg) * sc1a
  scfac(ircpg)   = sc1a


  !..n14 + p 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                              &
       &             zion(in14),aion(in14),zion(ih1),aion(ih1),    &
       &             jscr,screen_init,sc1a)

  ratdum(irnpg)  = ratraw(irnpg) * sc1a
  scfac(irnpg)   = sc1a


  !..o16 + p reactions
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(io16),aion(io16),zion(ih1),aion(ih1),   &
       &             jscr,screen_init,sc1a)

  ratdum(iropg)  = ratraw(iropg)
  scfac(iropg)   = 1.0e0


  !..the always fun triple alpha
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ihe4),aion(ihe4),4.0e0,8.0e0,           &
       &             jscr,screen_init,sc2a)

  sc3a          = sc1a * sc2a         

  ratdum(ir3a)   = ratraw(ir3a) * sc3a
  scfac(ir3a)    = sc3a

  ratdum(irg3a)  = ratraw(irg3a)
  scfac(irg3a)   = 1.0e0


  !..n14 + alpha 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(in14),aion(in14),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irnag)  = ratraw(irnag) * sc1a
  scfac(irnag)   = sc1a


  !..c12 + alpha 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(ircag)  = ratraw(ircag) * sc1a
  scfac(ircag)   = sc1a

  ratdum(iroga)  = ratraw(iroga)
  scfac(iroga)   = 1.0e0


  !..c12 + c12
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ic12),aion(ic12),zion(ic12),aion(ic12), &
       &             jscr,screen_init,sc1a)

  ratdum(ir1212) = ratraw(ir1212) * sc1a
  scfac(ir1212)  = sc1a


  !..c12 + o16
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ic12),aion(ic12),zion(io16),aion(io16), &
       &             jscr,screen_init,sc1a)

  ratdum(ir1216) = ratraw(ir1216) * sc1a
  scfac(ir1216)  = sc1a


  !..c12 + o16
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(io16),aion(io16),zion(io16),aion(io16), &
       &             jscr,screen_init,sc1a)

  ratdum(ir1616) = ratraw(ir1616) * sc1a
  scfac(ir1616)  = sc1a


  !..o16 + alpha 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(iroag)  = ratraw(iroag) * sc1a
  scfac(iroag)   = sc1a

  ratdum(irnega) = ratraw(irnega)
  scfac(irnega)  = 1.0e0


  !..ne20 + alpha 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(ine20),aion(ine20),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irneag) = ratraw(irneag) * sc1a
  scfac(irneag)  = sc1a

  ratdum(irmgga) = ratraw(irmgga)
  scfac(irmgga)  = 1.0e0


  !..mg24 to si28 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(img24),aion(img24),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irmgag) = ratraw(irmgag) * sc1a
  scfac(irmgag)  = sc1a

  ratdum(irsiga) = ratraw(irsiga)
  scfac(irsiga)  = 1.0e0

  ratdum(irmgap) = ratraw(irmgap) * sc1a
  scfac(irmgap)  = sc1a


  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             13.0e0,27.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(iralpa) = ratraw(iralpa) * sc1a
  scfac(iralpa)  = sc1a

  ratdum(iralpg) = ratraw(iralpg) * sc1a
  scfac(iralpg)  = sc1a

  ratdum(irsigp) = ratraw(irsigp)
  scfac(irsigp)  = 1.0e0


  !..si28 to s32
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(isi28),aion(isi28),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irsiag) = ratraw(irsiag) * sc1a
  scfac(irsiag)  = sc1a

  ratdum(irsga)  = ratraw(irsga)
  scfac(irsga)   = 1.0e0

  ratdum(irsiap) = ratraw(irsiap) * sc1a
  scfac(irsiap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             15.0e0,31.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irppa)  = ratraw(irppa) * sc1a
  scfac(irppa)   = sc1a

  ratdum(irppg)  = ratraw(irppg) * sc1a
  scfac(irppg)   = sc1a 

  ratdum(irsgp)  = ratraw(irsgp)
  scfac(irsgp)   = 1.0e0


  !..s32 to ar36
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(is32),aion(is32),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irsag)  = ratraw(irsag) * sc1a
  scfac(irsag)   = sc1a

  ratdum(irarga) = ratraw(irarga)
  scfac(irarga)  = 1.0e0

  ratdum(irsap)  = ratraw(irsap) * sc1a
  scfac(irsap)   = sc1a


  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             17.0e0,35.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irclpa) = ratraw(irclpa) * sc1a
  scfac(irclpa)  = sc1a 

  ratdum(irclpg) = ratraw(irclpg) * sc1a
  scfac(irclpg)  = sc1a

  ratdum(irargp) = ratraw(irargp)
  scfac(irargp)  = 1.0e0


  !..ar36 to ca40
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(iar36),aion(iar36),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irarag) = ratraw(irarag) * sc1a
  scfac(irarag)  = sc1a

  ratdum(ircaga) = ratraw(ircaga)
  scfac(ircaga)  = 1.0e0

  ratdum(irarap) = ratraw(irarap) * sc1a
  scfac(irarap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             19.0e0,39.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irkpa)  = ratraw(irkpa) * sc1a
  scfac(irkpa)   = sc1a

  ratdum(irkpg)  = ratraw(irkpg) * sc1a
  scfac(irkpg)   = sc1a

  ratdum(ircagp) = ratraw(ircagp)
  scfac(ircagp)  = 1.0e0


  !..ca40 to ti44
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(ircaag) = ratraw(ircaag) * sc1a
  scfac(ircaag)  = sc1a

  ratdum(irtiga) = ratraw(irtiga)
  scfac(irtiga)  = 1.0e0

  ratdum(ircaap) = ratraw(ircaap) * sc1a
  scfac(ircaap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             21.0e0,43.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irscpa) = ratraw(irscpa) * sc1a
  scfac(irscpa)  = sc1a

  ratdum(irscpg) = ratraw(irscpg) * sc1a
  scfac(irscpg)  = sc1a

  ratdum(irtigp) = ratraw(irtigp)
  scfac(irtigp)  = 1.0e0


  !..ti44 to cr48
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irtiag) = ratraw(irtiag) * sc1a
  scfac(irtiag)  = sc1a

  ratdum(ircrga) = ratraw(ircrga)
  scfac(ircrga)  = 1.0e0

  ratdum(irtiap) = ratraw(irtiap) * sc1a
  scfac(irtiap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             23.0e0,47.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irvpa)  = ratraw(irvpa) * sc1a
  scfac(irvpa)   = sc1a

  ratdum(irvpg)  = ratraw(irvpg) * sc1a
  scfac(irvpg)   = sc1a

  ratdum(ircrgp) = ratraw(ircrgp)
  scfac(ircrgp)  = 1.0e0


  !..cr48 to fe52
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(icr48),aion(icr48),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(ircrag) = ratraw(ircrag) * sc1a
  scfac(ircrag)  = sc1a

  ratdum(irfega) = ratraw(irfega)
  scfac(irfega)  = 1.0e0

  ratdum(ircrap) = ratraw(ircrap) * sc1a
  scfac(ircrap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             25.0e0,51.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(irmnpa) = ratraw(irmnpa) * sc1a
  scfac(irmnpa)  = sc1a

  ratdum(irmnpg) = ratraw(irmnpg) * sc1a
  scfac(irmnpg)  = sc1a

  ratdum(irfegp) = ratraw(irfegp)
  scfac(irfegp)  = 1.0e0


  !..fe52 to ni56
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(ife52),aion(ife52),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(irfeag) = ratraw(irfeag) * sc1a
  scfac(irfeag)  = sc1a

  ratdum(irniga) = ratraw(irniga)
  scfac(irniga)  = 1.0e0

  ratdum(irfeap) = ratraw(irfeap) * sc1a
  scfac(irfeap)  = sc1a

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             27.0e0,55.0e0,zion(ih1),aion(ih1),           &
       &             jscr,screen_init,sc1a)

  ratdum(ircopa) = ratraw(ircopa) * sc1a
  scfac(ircopa)  = sc1a

  ratdum(ircopg) = ratraw(ircopg) * sc1a
  scfac(ircopg)  = sc1a

  ratdum(irnigp) = ratraw(irnigp)
  scfac(irnigp)  = 1.0e0


  !..fe54
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ife52),aion(ife52),zion(ih1),aion(ih1), &
       &             jscr,screen_init,sc1a)

  ratdum(irfepg) = ratraw(irfepg)
  scfac(irfepg)  = 1.0e0

  ratdum(ircogp) = ratraw(ircogp)
  scfac(ircogp)  = 1.0e0



  ratdum(ir52ng) = ratraw(ir52ng)
  scfac(ir52ng)  = 1.0e0

  ratdum(ir53gn) = ratraw(ir53gn)
  scfac(ir53gn)  = 1.0e0

  ratdum(ir53ng) = ratraw(ir53ng)
  scfac(ir53ng)  = 1.0e0

  ratdum(ir54gn) = ratraw(ir54gn)
  scfac(ir54gn)  = 1.0e0

  ratdum(irheng) = ratraw(irheng)
  scfac(irheng)  = 1.0e0

  ratdum(irhegn) = ratraw(irhegn)
  scfac(irhegn)  = 1.0e0

  ratdum(irhng)  = ratraw(irhng)
  scfac(irhng)   = 1.0e0

  ratdum(irdgn)  = ratraw(irdgn)
  scfac(irdgn)   = 1.0e0


  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             1.0e0,2.0e0,zion(ih1),aion(ih1),             &
       &             jscr,screen_init,sc1a)

  ratdum(irdpg)  = ratraw(irdpg) * sc1a
  scfac(irdpg)   = sc1a


  ratdum(irhegp) = ratraw(irhegp)
  scfac(irhegp)  = 1.0e0

  ratdum(irpen)  = ratraw(irpen)
  scfac(irpen)   = 1.0e0

  ratdum(ispen)  = ratraw(ispen)
  scfac(ispen)   = 1.0e0

  ratdum(irnep)  = ratraw(irnep)
  scfac(irnep)   = 1.0e0

  ratdum(isnep)  = ratraw(isnep)
  scfac(isnep)   = 1.0e0

  ratdum(irn56ec) = ratraw(irn56ec)
  scfac(irn56ec)  = 1.0e0

  ratdum(isn56ec) = ratraw(isn56ec)
  scfac(isn56ec)  = 1.0e0

  ratdum(ifa)    = ratraw(ifa)
  scfac(ifa)     = 1.0e0

  ratdum(ifg)    = ratraw(ifg)
  scfac(ifg)     = 1.0e0


  !..reset the screen initialization flag
  screen_init = 0

  return
end subroutine bn_networkScreen

