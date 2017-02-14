!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Burn_data
!!
!! NAME
!!  Burn_data
!!
!! SYNOPSIS
!!
!!  use Burn_data
!!
!! DESCRIPTION 
!!  Burn_data is a fortran module that holds variables with
!!  the Burn Unit scope.  All variables located in Burn_data are
!!  accessible to subroutines in the Burn unit only.  (Note this 
!!  is a convention, there is nothing in the fortran language that
!!  prevents another routine from using Burn_data.  It is the FLASH3
!!  architecture that tries to enforce these rules.)
!!  
!!  All variables located in the Burn_data fortran module start with 
!! "bn_".  This is to indicate where they come from, and to make it easier
!! on the developer and user to see which variables are local variables and 
!! which belong to the Burn unit.
!!
!!
!! NOTES
!!  For those familiar with FLASH2:  The fortran data modules in FLASH3 
!!  (like this Burn_data) take the place of the database in FLASH2.
!!  In FLASH2 all the data was centrally located in the database.  In 
!!  FLASH3 each unit stores its own data.  When another unit needs to 
!!  get data from a different unit, accessor functions are used.  For example
!!  simulation time is stored in the Driver owned variable dr_simTime.  The
!!  IO unit queries the Driver unit with the accessor function
!!  Driver_getSimTime() to get the simTime and to check to see if a checkpoint
!!  needs to be written.
!! 
!! This file Burn_data was originally primarily network_common.fh in Flash2
!! Some common blocks have been moved out to other files
!!
!!***

Module Burn_data

  ! This variable changes with the network size and therefore must be initialized in each aprox directory
  use bn_dataNetworkSize, ONLY:  nrat, nratp1

  implicit none

#include "constants.h"
#include "Flash.h"

! Runtime Parameters
  integer, save :: bn_algebra, bn_odeStepper
  integer, save :: bn_meshMe

  logical, save :: bn_useBurn  
  logical, save :: bn_useShockBurn, bn_useBurnTable
  real, save    :: bn_smallx, bn_enucDtFactor
  real, save    :: bn_nuclearTempMin, bn_nuclearTempMax, bn_nuclearDensMin, &
       &           bn_nuclearDensMax, bn_nuclearNI56Max

! were in network_size.fh, specific to aprox13, I believe.
!! Now in bn_dataNetworkSize
!..this sets the size of the number of nuclear reaction rates
!..only critical to the burn modules
!..nrat    = number of reaction rates
!!  integer, parameter :: nrat = 59
!!  integer, parameter :: nratp1 = nrat+1


!..nuclear reation network common block declarations

!..xmass   = mass fractions
!..ymass   = molar fractions
!..aion    = number of nucleons
!..aioninv  = 1/aion
!..zion    = number of protons
!..zionsq  = zion*zion
!..bion    = binding energies
!..ionnam  = name of isotope

!..ratnam  = name of reaction rates
!..ratdum  = the raw reaction rates (unscreened) as an array
!..ratdum  = the screened reaction rates as an array

!..ra1     = nucleons in one reacting channel
!..rz1     = charge in one reacting channel
!..ra2     = nucleons in other reacting channel
!..rz2     = charge in other reacting channel
!..zs13    = (z1+z2)**(1./3.)
!..zhat    = combination of z1 and z2 raised to the 5/3 power
!..zhat2   = combination of z1 and z2 raised to the 5/12 power
!..lzav    = log of effective charge
!..aznut   = combination of a1,z1,a2,z2 raised to 1/3 power
!..scfac   = screening factors

!..sneut   = total neutrino energy losss rate
!..sphot   = neutrino losss from photodisintegration
!..splas   = neutrino loss from plasmons
!..spair   = neutrino loss from pair production
!..sbrem   = neutrino less from bremmstrahlung
!..srecomb = neutrino loss from recombination

!..xoktot  = total number of burner steps taken
!..xbadtot = total number of bad, but redone, burner steps
!..xkbrn   = total number of times the burner was called

! were in common block /netc2/
  character(len=5), save ::    ratnam(nrat)
  character(len=4), save ::    ionam(NSPECIES)
  
! were in common block /netc4/
  real,save,dimension(NSPECIES) :: xmass,ymass,aion,zion,bion,aioninv,zionsq    ! used to have aionin too
  real,save        ::                                                   &
     &                 rz1(nrat+1),ra1(nrat+1),                         &
     &                 rz2(nrat+1),ra2(nrat+1),                         &
     &                 zs13(nrat),zhat(nrat),zhat2(nrat),               &
     &                 lzav(nrat),aznut(nrat),scfac(nrat),              &
     &                 zs13inv(nrat),                                   &
     &                 ratraw(nrat),ratdum(nrat),                       &
     &                 xoktot,xbadtot,xkbrn                           
      integer,save      ::  isflag(nrat+1)

! were in common block /netc5/
      real, save ::    sneut,sphot,spair,splas,sbrem,srecomb



!..for tabular evaluation of the raw reaction rates
!..allow storage for 120/points per decade
!..logical useBurnTables for determing if tables are to be used
! were in common block /rcm2t/
      integer, parameter   ::        nrattab = 481
      real, save      ::        rattab(nrat,nrattab),                   & 
           ttab(nrattab),dtab(nrat)


!..for nice identification of 62 key isotopes
!were in common block /netc8/

       integer, parameter ::   nisotp=64
       integer, save ::                                                 &
     &         ih1,iprot,ineut,ihe4,ih2,ih3,ihe3,ili6,ili7,ibe7,ibe9,   &
     &         ib8,ib10,ib11,ic11,ic12,ic13,ic14,in13,in14,in15,io14,   &
     &         io15,io16,io17,io18,if17,if18,if19,ine18,ine19,ine20,    &
     &         ine21,ine22,ina21,ina22,ina23,img22,img23,img24,img25,   &
     &         img26,ial25,ial26,ial27,isi27,isi28,isi29,isi30,ip30,    &
     &         ip31,is30,is31,is32,iar36,ica40,iti44,icr48,ife52,       &
     &         ife54,ini56,izn60,ifuel,iash
     !!  This equivalent array is defined for ease of initialization, etc. 
     !!  Used in bn_initNetwork
     integer, dimension(nisotp), save ::    isotp
     equivalence      (isotp(1),ih1)


!..for easy aprox13 rate identification:
! were in common block /netc12/, now in file Aprox13_data

!..add these rates for the aprox19 network (but appear to be also used in Aprox13)
! were in common block /netc13/, now in file Aprox13_data

!..add these rates for the pp123 network
! were in common block /netc13a/, now in file Burn_dataPP123

!..add these rates for the cno network
! were in common block /netc13b/, now in file Burn_dataCNO

!..add these rates for the rp network
! were in common block /netc13c/, now file Burn_dataRP

!..add these rates for the generic 2 fluid network
! was in common block /netc13d/, didn't add to separate file
      integer, save ::          irfuel   

end Module Burn_data
