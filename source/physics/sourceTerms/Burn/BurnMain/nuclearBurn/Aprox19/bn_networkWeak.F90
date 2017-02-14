!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/Aprox19/bn_networkWeak
!!
!! NAME
!!
!!  bn_network
!!
!!
!! SYNOPSIS
!!
!!  call bn_networkWeak(real(IN) ::y)
!!
!!
!! DESCRIPTION
!!
!!  aprox19 was the third network installed in flash 1.5
!!  it is like the original iso13 network, but adds protons and he3
!!  for pp burning, nitrogen to do cno and hot cno (beta limited)
!!  burning, fe54 and photodisintegration protons and neutrons
!!
!! ARGUMENTS
!!
!!  y(NSPECIES) -- mass of the NSPECIES
!!
!! NOTES
!!
!!  uses
!!  bn_initNetwork to initialize the aprox19 network
!!  bn_burner to drive the aprox19 network
!!
!!
!!
!! routine aprox19 sets up the odes 
!! routine bn_networkDenseJakob sets up the dense aprox19 jacobian
!! routine bn_networkSparsePointers builds the nonzero locations for bn_networkSparseJakob
!! routine bn_networkSparseJakob sets up the sparse aprox19 jacobian 
!! routine bn_networkRates generates the raw reaction rates for routine aprox19
!! routine bn_networkTable generates the raw rates using table interpolation
!! routine bn_networkScreen applies screening corrections to the raw rates
!!***

subroutine bn_networkWeak(y)

#include "constants.h"   
#include "Flash.h"
#include "Eos.h"
#include "Eos_components.h"

  use Eos_interface, ONLY: Eos
  use bn_interface, ONLY: bn_ecapnuc, bn_mazurek

  use Burn_dataEOS, ONLY: btemp, den, ye
  use Burn_data
  use bn_dataAprox19

  implicit none

  ! it appears that den and btemp are passed in via eos_common.fh or Burn_dataEOS
  ! we use the declarations of abar and zbar from there as well. 

  !..electron capture rates on nucleons for aprox19
  !..note they are composition dependent

  !..declare
  real, intent(IN) :: y(NSPECIES)

  integer   i, specieMap
  real      xn(NSPECIES),                    &
       &          dpt,dpd,det,ded,gammac,pel,xne,eta, &
       &          rpen,rnep,spen,snep, c_v, c_p

  !! Needed for Eos call
  integer :: vecLen
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  real, dimension(EOS_NUM) :: eosData
  real, dimension(NSPECIES) :: massFrac  ! in the order required by the rest of Flash!


  !..generate the mass fractions from the passed composition
  real :: entropy, dst, dsd
  do i=1,NSPECIES
     xn(i) = y(i)*aion(i)
  enddo

  !..get the degeneracy parameter eta
  vecLen = 1
  mask = .true.

  mask(EOS_DEA) = .FALSE.
  mask(EOS_DEZ) = .FALSE.
  
  eosData(EOS_DENS) = den
#ifdef FLASH_UHD_3T
#if EOSCOMP_MATTER == EOSCOMP_ELE
  eosData(EOS_TEMPELE) = btemp
  eosData(EOS_TEMPION) = 0.0
#else
  eosData(EOS_TEMPELE) = btemp
  eosData(EOS_TEMPION) = btemp
#endif
  eosData(EOS_TEMPRAD) = btemp
  eosData(EOS_TEMP) = btemp
#else
  eosData(EOS_TEMP) = btemp
#endif

  ! perform mapping from random burn order to Flash order
  do i = 1, NSPECIES
     call bn_mapNetworkToSpecies(i,specieMap)
     massFrac(specieMap - SPECIES_BEGIN + 1) = xn(i)
  end do

#ifdef FLASH_UHD_3T
        call Eos(MODE_DENS_TEMP_GATHER,vecLen,eosData,massFrac,mask)
#else
        call Eos(MODE_DENS_TEMP,vecLen,eosData,massFrac,mask)
#endif

  eta = eosData(EOS_ETA)


  !..get the electron capture rates
  call bn_ecapnuc(eta,btemp,rpen,rnep,spen,snep)
  ratraw(irpen)   = rpen
  ratraw(irnep)   = rnep
  ratraw(ispen)   = spen
  ratraw(isnep)   = snep


  !..ni56 electron capture rate
  call bn_mazurek(btemp,den,y(ini56),ye,                &
                    ratdum(irn56ec),ratdum(isn56ec))

  return
end subroutine bn_networkWeak



