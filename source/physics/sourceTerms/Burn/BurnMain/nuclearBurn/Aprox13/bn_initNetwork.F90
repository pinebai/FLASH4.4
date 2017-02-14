!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox13/bn_initNetwork
!!
!! NAME
!!  
!!  bn_initNetwork
!!
!!
!! SYNOPSIS
!! 
!!  bn_initNetwork()
!!
!!  
!! DESCRIPTION
!!
!!  Initializes variables and integers used in the network approximation Aprox13.
!!  Called from Burn_init.
!!
!!  this is an alpha chain + heavy ion network with (a,p)(p,g) links
!!
!!  isotopes: he4,  c12,  o16,  ne20, mg24, si28, s32,
!!            ar36, ca40, ti44, cr48, fe52, ni56
!!
!!***


subroutine bn_initNetwork

  use Burn_data
  use bn_dataAprox13

#include "Flash.h"
  implicit none

  !..this routine initializes stuff for the aprox13 network
  !  Nearly all variables defined in the data structures.

  integer ::   i


  !..zero all the isotope pointers
  do i=1,nisotp
     isotp(i)   = 0
  enddo

  !..zero the steps taken
  xoktot  = 0.0e0
  xbadtot = 0.0e0
  xkbrn   = 0.0e0


  !..set the id numbers of the elements
  !! These need to correspond to what is defined in the Flash.h!
  !! However, this whole thing doesn't seem to work independent of the hardwired numbers
  ihe4  = 1
  ic12  = 2
  io16  = 3
  ine20 = 4
  img24 = 5
  isi28 = 6
  is32  = 7
  iar36 = 8
  ica40 = 9
  iti44 = 10
  icr48 = 11
  ife52 = 12
  ini56 = 13 



  !..set the names of the elements
  ionam(ihe4)  = 'he4 '
  ionam(ic12)  = 'c12 '
  ionam(io16)  = 'o16 '
  ionam(ine20) = 'ne20'
  ionam(img24) = 'mg24'
  ionam(isi28) = 'si28'
  ionam(is32)  = 's32 '
  ionam(iar36) = 'ar36'
  ionam(ica40) = 'ca40'
  ionam(iti44) = 'ti44'
  ionam(icr48) = 'cr48'
  ionam(ife52) = 'fe52'
  ionam(ini56) = 'ni56'

  !..set the number of nucleons in the element
  aion(ihe4)  = 4.0e0
  aion(ic12)  = 12.0e0   
  aion(io16)  = 16.0e0   
  aion(ine20) = 20.0e0   
  aion(img24) = 24.0e0  
  aion(isi28) = 28.0e0  
  aion(is32)  = 32.0e0  
  aion(iar36) = 36.0e0  
  aion(ica40) = 40.0e0  
  aion(iti44) = 44.0e0  
  aion(icr48) = 48.0e0  
  aion(ife52) = 52.0e0  
  aion(ini56) = 56.0e0  
  do i=1,NSPECIES
     aioninv(i) = 1.0e0/aion(i)
  enddo


  !..set the number of protons in the element
  zion(ihe4)  = 2.0e0
  zion(ic12)  = 6.0e0
  zion(io16)  = 8.0e0
  zion(ine20) = 10.0e0   
  zion(img24) = 12.0e0  
  zion(isi28) = 14.0e0  
  zion(is32)  = 16.0e0  
  zion(iar36) = 18.0e0  
  zion(ica40) = 20.0e0  
  zion(iti44) = 22.0e0  
  zion(icr48) = 24.0e0  
  zion(ife52) = 26.0e0  
  zion(ini56) = 28.0e0  
  do i=1,NSPECIES
     zionsq(i) = zion(i) * zion(i)
  enddo


  !..set the binding energy of the element
  bion(ihe4)  =  28.29603e0 
  bion(ic12)  =  92.16294e0
  bion(io16)  = 127.62093e0 
  bion(ine20) = 160.64788e0 
  bion(img24) = 198.25790e0 
  bion(isi28) = 236.53790e0 
  bion(is32)  = 271.78250e0 
  bion(iar36) = 306.72020e0 
  bion(ica40) = 342.05680e0  
  bion(iti44) = 375.47720e0 
  bion(icr48) = 411.46900e0 
  bion(ife52) = 447.70800e0
  bion(ini56) = 484.00300e0 


  !..set the id numbers of the reaction rates
  ir3a   = 1
  irg3a  = 2
  ircag  = 3
  ir1212 = 4
  ir1216 = 5
  ir1616 = 6
  iroga  = 7
  iroag  = 8
  irnega = 9
  irneag = 10
  irmgga = 11
  irmgag = 12
  irsiga = 13
  irmgap = 14
  iralpa = 15
  iralpg = 16
  irsigp = 17
  irsiag = 18
  irsga  = 19
  irsiap = 20
  irppa  = 21
  irppg  = 22
  irsgp  = 23
  irsag  = 24
  irarga = 25
  irsap  = 26
  irclpa = 27
  irclpg = 28
  irargp = 29
  irarag = 30
  ircaga = 31
  irarap = 32
  irkpa  = 33
  irkpg  = 34
  ircagp = 35
  ircaag = 36
  irtiga = 37
  ircaap = 38
  irscpa = 39
  irscpg = 40
  irtigp = 41
  irtiag = 42
  ircrga = 43
  irtiap = 44
  irvpa  = 45
  irvpg  = 46
  ircrgp = 47
  ircrag = 48
  irfega = 49
  ircrap = 50
  irmnpa = 51
  irmnpg = 52
  irfegp = 53
  irfeag = 54
  irniga = 55
  irfeap = 56
  ircopa = 57
  ircopg = 58
  irnigp = 59

  !..set the names of the reaction rates
  ratnam(ir3a)   = 'r3a  '
  ratnam(irg3a)  = 'rg3a '
  ratnam(ircag)  = 'rcag '
  ratnam(ir1212) = 'r1212'
  ratnam(ir1216) = 'r1216'
  ratnam(ir1616) = 'r1616'
  ratnam(iroga)  = 'roga '
  ratnam(iroag)  = 'roag '
  ratnam(irnega) = 'rnega'
  ratnam(irneag) = 'rneag'
  ratnam(irmgga) = 'rmgga'
  ratnam(irmgag) = 'rmgag'
  ratnam(irsiga) = 'rsiga'
  ratnam(irmgap) = 'rmgap'
  ratnam(iralpa) = 'ralpa'
  ratnam(iralpg) = 'ralpg'
  ratnam(irsigp) = 'rsigp'
  ratnam(irsiag) = 'rsiag'
  ratnam(irsga)  = 'rsga '
  ratnam(irsiap) = 'rsiap'
  ratnam(irppa)  = 'rppa '
  ratnam(irppg)  = 'rppg '
  ratnam(irsgp)  = 'rsgp '
  ratnam(irsag)  = 'rsag '
  ratnam(irarga) = 'rarga'
  ratnam(irsap)  = 'rsap '
  ratnam(irclpa) = 'rclpa'
  ratnam(irclpg) = 'rclpg'
  ratnam(irargp) = 'rargp'
  ratnam(irarag) = 'rarag'
  ratnam(ircaga) = 'rcaga'
  ratnam(irarap) = 'rarap'
  ratnam(irkpa)  = 'rkpa '
  ratnam(irkpg)  = 'rkpg '
  ratnam(ircagp) = 'rcagp'
  ratnam(ircaag) = 'rcaag'
  ratnam(irtiga) = 'rtiga'
  ratnam(ircaap) = 'rcaap'
  ratnam(irscpa) = 'rscpa'
  ratnam(irscpg) = 'rscpg'
  ratnam(irtigp) = 'rtigp'
  ratnam(irtiag) = 'rtiag'
  ratnam(ircrga) = 'rcrga'
  ratnam(irtiap) = 'rtiap'
  ratnam(irvpa)  = 'rvpa '
  ratnam(irvpg)  = 'rvpg '
  ratnam(ircrgp) = 'rcrgp'
  ratnam(ircrag) = 'rcrag'
  ratnam(irfega) = 'rfega'
  ratnam(ircrap) = 'rcrap'
  ratnam(irmnpa) = 'rmnpa'
  ratnam(irmnpg) = 'rmnpg'
  ratnam(irfegp) = 'rfegp'
  ratnam(irfeag) = 'rfeag'
  ratnam(irniga) = 'rniga'
  ratnam(irfeap) = 'rfeap'
  ratnam(ircopa) = 'rcopa'
  ratnam(ircopg) = 'rcopg'
  ratnam(irnigp) = 'rnigp'

  return
end subroutine bn_initNetwork
