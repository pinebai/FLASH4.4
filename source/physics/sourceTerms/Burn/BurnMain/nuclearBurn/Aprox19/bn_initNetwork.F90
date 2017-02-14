!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_initNetwork
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
!!  aprox19 was the third network installed in flash 1.5
!!  it is like the original iso13 network, but adds protons and he3
!!  for pp burning, nitrogen to do cno and hot cno (beta limited)
!!  burning, fe54 and photodisintegration protons and neutrons
!!
!!***

subroutine bn_initNetwork

  use Burn_data
  use bn_dataAprox19

  implicit none

#include "Flash.h"

  !..
  !..this routine initializes stuff for the aprox19 network


  !..declare local variables
   integer :: i
!  integer          i,isotp(nisotp)
!  equivalence      (isotp(1),ih1)



  !..zero all the isotope pointers
  real :: entropy, dst, dsd
  do i=1,nisotp
     isotp(i)   = 0
  enddo


  !..zero the steps taken
  xoktot  = 0.0e0
  xbadtot = 0.0e0
  xkbrn   = 0.0e0


  !..set the id numbers of the elements
  ih1   = 1
  ihe3  = 2
  ihe4  = 3
  ic12  = 4
  in14  = 5
  io16  = 6
  ine20 = 7
  img24 = 8
  isi28 = 9
  is32  = 10
  iar36 = 11
  ica40 = 12
  iti44 = 13
  icr48 = 14
  ife52 = 15
  ife54 = 16
  ini56 = 17 
  ineut = 18
  iprot = 19


  !..set the names of the elements
  ionam(ih1)   = 'h1  '
  ionam(ihe3)  = 'he3 '
  ionam(ihe4)  = 'he4 '
  ionam(ic12)  = 'c12 '
  ionam(in14)  = 'n14 '
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
  ionam(ife54) = 'fe54'
  ionam(ini56) = 'ni56'
  ionam(ineut) = 'neut'
  ionam(iprot) = 'prot'

  !..set the number of nucleons in the element
  aion(ih1)   = 1.0e0
  aion(ihe3)  = 3.0e0
  aion(ihe4)  = 4.0e0
  aion(ic12)  = 12.0e0   
  aion(in14)  = 14.0e0   
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
  aion(ife54) = 54.0e0  
  aion(ini56) = 56.0e0  
  aion(ineut) = 1.0e0
  aion(iprot) = 1.0e0
  do i=1,NSPECIES
     aioninv(i) = 1.0e0/aion(i)
  enddo


  !..set the number of protons in the element
  zion(ih1)   = 1.0e0
  zion(ihe3)  = 2.0e0
  zion(ihe4)  = 2.0e0
  zion(ic12)  = 6.0e0
  zion(in14)  = 7.0e0
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
  zion(ife54) = 26.0e0  
  zion(ini56) = 28.0e0  
  zion(ineut) = 0.0e0
  zion(iprot) = 1.0e0
  do i=1,NSPECIES
     zionsq(i) = zion(i) * zion(i)
  enddo

  !..set the binding energy of the element
  bion(ih1)   = 0.0e0
  bion(ihe3)  = 7.71819e0
  bion(ihe4)  = 28.29603e0 
  bion(ic12)  = 92.16294e0
  bion(in14)  = 104.65998e0
  bion(io16)  = 127.62093e0 
  bion(ine20) = 160.64788e0 
  bion(img24) = 198.2579e0 
  bion(isi28) = 236.5379e0 
  bion(is32)  = 271.7825e0 
  bion(iar36) = 306.7202e0 
  bion(ica40) = 342.0568e0  
  bion(iti44) = 375.4772e0 
  bion(icr48) = 411.469e0 
  bion(ife52) = 447.708e0
  bion(ife54) = 471.7696e0
  bion(ini56) = 484.003e0 
  bion(ineut) = 0.0e0
  bion(iprot) = 0.0e0

  !..set the id numbers of the reaction rates
  irpp   = 1
  ir33   = 2
  ir34   = 3
  ircpg  = 4
  irnpg  = 5
  iropg  = 6
  ir3a   = 7
  irg3a  = 8
  irnag  = 9
  ircag  = 10
  ir1212 = 11
  ir1216 = 12
  ir1616 = 13
  iroga  = 14
  iroag  = 15
  irnega = 16
  irneag = 17
  irmgga = 18
  irmgag = 19
  irsiga = 20
  irmgap = 21
  iralpa = 22
  iralpg = 23
  irsigp = 24
  irsiag = 25
  irsga  = 26
  irsiap = 27
  irppa  = 28
  irppg  = 29
  irsgp  = 30
  irsag  = 31
  irarga = 32
  irsap  = 33
  irclpa = 34
  irclpg = 35
  irargp = 36
  irarag = 37
  ircaga = 38
  irarap = 39
  irkpa  = 40
  irkpg  = 41
  ircagp = 42
  ircaag = 43
  irtiga = 44
  ircaap = 45
  irscpa = 46
  irscpg = 47
  irtigp = 48
  irtiag = 49
  ircrga = 50
  irtiap = 51
  irvpa  = 52
  irvpg  = 53
  ircrgp = 54
  ircrag = 55
  irfega = 56
  ircrap = 57
  irmnpa = 58
  irmnpg = 59
  irfegp = 60
  irfeag = 61
  irniga = 62
  irfeap = 63
  ircopa = 64
  ircopg = 65
  irnigp = 66
  irfepg = 67
  ircogp = 68
  ir52ng = 69
  ir53gn = 70
  ir53ng = 71
  ir54gn = 72
  irheng = 73
  irhegn = 74
  irhng  = 75
  irdgn  = 76
  irdpg  = 77
  irhegp = 78

  irpen   = 79
  ispen   = 80
  irnep   = 81
  isnep   = 82
  irn56ec = 83
  isn56ec = 84
  ifa     = 85
  ifg     = 86


  !..set the names of the reaction rates
  ratnam(irpp)   = 'rpp  '
  ratnam(ir33)   = 'r33  '
  ratnam(ir34)   = 'r34  '
  ratnam(ircpg)  = 'rcpg '
  ratnam(irnpg)  = 'rnpg '
  ratnam(iropg)  = 'ropg '
  ratnam(ir3a)   = 'r3a  '
  ratnam(irg3a)  = 'rg3a '
  ratnam(irnag)  = 'rnag '
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
  ratnam(irfepg) = 'rfepg'
  ratnam(ircogp) = 'rcogp'
  ratnam(ir52ng) = 'r52ng'
  ratnam(ir53gn) = 'r53gn'
  ratnam(ir53ng) = 'r53ng'
  ratnam(ir54gn) = 'r54gn'
  ratnam(irheng) = 'rheng'
  ratnam(irhegn) = 'rhegn'
  ratnam(irhng)  = 'rhng '
  ratnam(irdgn)  = 'rdgn '
  ratnam(irdpg)  = 'rdpg '
  ratnam(irhegp) = 'rhegp'

  ratnam(irpen)   = 'rpen '
  ratnam(ispen)   = 'spen '
  ratnam(irnep)   = 'rnep '
  ratnam(isnep)   = 'snep '
  ratnam(irn56ec) = 'r56ec'
  ratnam(isn56ec) = 's56ec'
  ratnam(ifa)     = 'rfa  '
  ratnam(ifg)     = 'rfg  '


  return 
end subroutine bn_initNetwork

