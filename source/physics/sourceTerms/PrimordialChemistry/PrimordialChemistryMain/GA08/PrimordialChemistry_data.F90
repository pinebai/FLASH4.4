!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/PrimordialChemistry_data
!!
!! NAME
!!
!!  PrimordialChemistry_data
!!
!!
!! SYNOPSIS
!!
!! use PrimordialChemistry_data
!!
!! DESCRIPTION
!!
!! Store the data for the Globular PrimordialChemistry unit
!!
!! ARGUMENTS
!!***

module PrimordialChemistry_data
    
  use pchem_dataNetworkSize, ONLY: nrat, nratp1

  implicit none
  !! Put all the required data here
      
#include "constants.h"
#include "Flash.h"
! There are probably some runtime stuff that should be here. 
! since I am following Burn and they have some. But for now
! I am going to put this off till later.
! Below will be a lot of stuff from Burn
!RUNTIME PARAMETERS NEED TO GO HERE, BUT DON'T KNOW WHAT!!

  !Runtime Paramters
  integer, save :: pchem_algebra, pchem_odeStepper
  integer, save :: pchem_meshMe

  logical, save :: pchem_usePrimordialChemistry
  logical, save :: pchem_useShockBurn




!..xmass   = mass fractions
!..ymass   = molar fractions
!..aion    = number of nucleons
!..aioninv = 1/aion
!..zion    = number of protons
!..zionsq  = zion*zion
!..bion    = binding energies
!..ionnam  = name of ion

!..ratnam  = name of reaction rate
!..ratdum  = the raw reaction rates (unscreened) as an array
!..ratdum  = the screened reaction rates as an array

!..ra1     = nucleons in one reacting channel
!..rz1     = charge in one reacting channel
!..ra2     = nucleons in other reacting channel
!..rz2     = charge in other reacting channel

character(len=11), save ::      ratnam(nrat)
character(len=4) , save ::      ionam(NSPECIES)

real,save,dimension(NSPECIES)  :: xmass,ymass,aion,zion,bion,aioninv,zionsq
real,save                      :: rz1(nrat+1), ra1(nrat+1), rz2(nrat+1), ra2(nrat+1), &
                                  zs13(nrat), zhat(nrat), zhat2(nrat), lzav(nrat), aznut(nrat), &
                                  scfac(nrat), zs13inv(nrat),ratraw(nrat), ratdum(nrat), xoktot, &
                                  xbadtot, xkburn, zz(111)
integer,save            :: isflag(nrat+1)

integer,parameter ::    nrattab = 481  !Not sure what this will do, need to ask
real,save        ::     rattab(nrat,nrattab), ttab(nrattab), dtab(nrat)

!!..Define all the ions as integers
integer, parameter :: nisotp = 16
integer,save :: iHP, iH, iHM, iDP, iD, iDM, iHEPP, iHEP, iHE, iH2P, iH2, iHDP, iHD, iELEC, iD2, iD2P
integer, dimension(nisotp), save :: isotp
equivalence (isotp(1),ihp)
 integer, save :: irfuel  !not sure if I will need this

real, save  :: pchem_fracHydrogen, pchem_fracHelium, pchem_fracDeuterium
real, save  :: pchem_j21, pchem_fshh2, pchem_fshhd
integer, save :: pchem_doCool, pchem_mCool, pchem_ccCase, pchem_rcCase


real, save :: amu
real, save :: pchem_tradmin, pchem_tradmax, pchem_dradmin, pchem_dradmax, pchem_massFracH, pchem_noCool
real, save, dimension(353) :: metal_cooling, hhe_cooling, will_hhe_cooling

end module PrimordialChemistry_data
