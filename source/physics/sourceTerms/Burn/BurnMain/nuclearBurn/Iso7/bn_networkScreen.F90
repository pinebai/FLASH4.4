!!****if* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_networkScreen
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
!!
!!***
subroutine bn_networkScreen(y)
  
  use bn_interface,ONLY: bn_screen4

  use Burn_dataEOS, ONLY: abar, zbar, z2bar, ytot1
  use Burn_data
  use bn_dataIso7

#include "Flash.h"

  implicit none



  !..this routine computes the screening factors
  !..and applies them to the raw reaction rates,
  !..producing the final reaction rates used by the
  !..right hand sides and jacobian matrix elements

  !..declare
  real, intent(IN) :: y(1)
  !! local
  integer          i,j,k,jscr,screen_init,screen_on
  parameter        (screen_on = 1)
  real             sc1a,sc2a,sc3a,zbarxx,z2barxx

  !..initialize
  data             screen_init/1/


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


  !..first the always fun triple alpha and its inverse
  jscr = 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ihe4),aion(ihe4),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ihe4),aion(ihe4),4.0e0,8.0e0,           &
       &             jscr,screen_init,sc2a)

  sc3a          = sc1a * sc2a                            

  ratdum(ir3a)  = ratraw(ir3a) * sc3a
  scfac(ir3a)   = sc3a

  ratdum(irg3a)  = ratraw(irg3a) 
  scfac(irg3a)   = 1.0e0


  !..c12 to o16 
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(ic12),aion(ic12),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(ircag)  = ratraw(ircag)  * sc1a
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


  !..o16 + o16
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(io16),aion(io16),zion(io16),aion(io16), &
       &             jscr,screen_init,sc1a)

  ratdum(ir1616) = ratraw(ir1616) * sc1a
  scfac(ir1616)  = sc1a


  !..o16 to ne20
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                             &
       &             zion(io16),aion(io16),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(iroag)  = ratraw(iroag) * sc1a 
  scfac(iroag)   = sc1a

  ratdum(irnega) = ratraw(irnega) 
  scfac(irnega)  = 1.0e0


  !..o16 to mg24
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
#if 0


  !..ca40 to ti44
  jscr = jscr + 1
  call bn_screen4(zbar,abar,z2bar,                               &
       &             zion(ica40),aion(ica40),zion(ihe4),aion(ihe4), &
       &             jscr,screen_init,sc1a)

  ratdum(ircaag) = ratraw(ircaag) * sc1a
  scfac(ircaag)  = sc1a

  ratdum(irtiga) = ratraw(irtiga) 
  scfac(irtiga)  = 1.0e0
#endif

  !..reset the screen initialization flag
  screen_init = 0

  return
end subroutine bn_networkScreen


