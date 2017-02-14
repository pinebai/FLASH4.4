!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Iso7/bn_initNetwork
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
!!  Initialize the Iso7 network, called from Burn_init.
!!     
!!  Iso7, the second network installed in flash 1.5
!!  like the original iso13 network, but adds
!!  (a,p)(p,g) sequences through fe52 to ni56 and neutrino losses
!!  Routine bn_burner drives the aprox13 network
!!     given time step tstep, temperature temp, density density, and
!!     composition xIn, this routine returns the burned composition xOut
!!     and the energy generation rate sdotRate.
!!
!!***
subroutine bn_initNetwork

  use Burn_data
  use bn_dataIso7

#include "Flash.h"

  implicit none

  !..
  !..this routine initializes stuff for the iso7 network
  !..
  !..declare
  integer :: i

!! Can't have an equivalence within a module
!!  integer          isotp(nisotp)
!!  equivalence      (isotp(1),ih1)



  !..zero all the isotope pointers
  do i=1,nisotp
     isotp(i)   = 0
  enddo


  !..zero the steps taken
  xoktot  = 0.0e0
  xbadtot = 0.0e0
  xkbrn   = 0.0e0



  !..set the id numbers of the elements
  ihe4  = 1
  ic12  = 2
  io16  = 3
  ine20 = 4
  img24 = 5
  isi28 = 6
  ini56 = 7 


  !..set the names of the elements
  ionam(ihe4)  = 'he4 '
  ionam(ic12)  = 'c12 '
  ionam(io16)  = 'o16 '
  ionam(ine20) = 'ne20'
  ionam(img24) = 'mg24'
  ionam(isi28) = 'si28'
  ionam(ini56) = 'ni56'


  !..set the number of nucleons in the element
  aion(ihe4)  = 4.0e0
  aion(ic12)  = 12.0e0   
  aion(io16)  = 16.0e0   
  aion(ine20) = 20.0e0   
  aion(img24) = 24.0e0  
  aion(isi28) = 28.0e0  
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
  bion(ini56) = 484.00300e0 



  !..set the id numbers of the reaction rates
  ircag  = 1
  iroga  = 2
  ir3a   = 3
  irg3a  = 4
  ir1212 = 5
  ir1216 = 6
  ir1616 = 7
  iroag  = 8
  irnega = 9
  irneag = 10
  irmgga = 11
  irmgag = 12
  irsiga = 13
  ircaag = 14
  irtiga = 15

  !..set the names of the reaction rates
  ratnam(ircag)  = 'rcag '
  ratnam(iroga)  = 'roga '
  ratnam(ir3a)   = 'r3a  '
  ratnam(irg3a)  = 'rg3a '
  ratnam(ir1212) = 'r1212'
  ratnam(ir1216) = 'r1216'
  ratnam(ir1616) = 'r1616'
  ratnam(iroag)  = 'roag '
  ratnam(irnega) = 'rnega'
  ratnam(irneag) = 'rneag'
  ratnam(irmgga) = 'rmgga'
  ratnam(irmgag) = 'rmgag'
  ratnam(irsiga) = 'rsiga'
  ratnam(ircaag) = 'rcaag'
  ratnam(irtiga) = 'rtiga'

  return
end subroutine bn_initNetwork
