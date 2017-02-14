!!****if* source/physics/Eos/EosNuclear/eos_nuclear
!!
!! NAME
!!
!!  eos_nuclear
!!
!! SYNOPSIS
!!  
!!  call eos_nuclear()
!!                 
!!
!! DESCRIPTION
!!
!!  Main driver routine for the nuclear EOS.  Calls the kernel 
!!  routine nuc_eos_short. 
!!
!! ARGUMENTS
!!
!! NOTES
!!      Parts of this unit are released under a different license than the
!!      usual FLASH license.  Specifically, some subroutines in the kernel 
!!      directory are released under the Creative Commons 
!!      attribution-noncommercial-share alike license.  Basically, if you use this
!!      unit in your work, the license requires that you cite the two articles 
!!      mentioned below.  More details may be found here:  
!!      stellarcollapse.org/equationofstate.
!!
!!      * O'Connor, E.P., & Ott, C.D. 2010, CQGra, 27, 114103
!!      * Couch, S.M. 2013, ApJ, 765, 29
!!
!!
!!***
subroutine eos_nuclear(mode,vecLen,eosData,massFrac,mask)

  use Driver_interface, ONLY : Driver_abortFlash
  use Multispecies_interface, ONLY : Multispecies_getSumInv, &
       Multispecies_getSumFrac
  use Logfile_interface, ONLY:  Logfile_stampMessage
  use Eos_data, ONLY : eos_singleSpeciesA, eos_singleSpeciesZ, eos_smalle, eos_smallt
  use eosmodule, ONLY : precision, temp_mev_to_kelvin, &
       eos_rhomin, eos_yemax, eos_tempmin, e_zeroPoint

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#ifdef FLASH_MULTISPECIES
#include "Multispecies.h"
#endif

  !     Arguments
  integer, INTENT(in) :: mode, vecLen
  real, INTENT(inout), dimension(vecLen*EOS_NUM) :: eosData
  real, optional,INTENT(in), dimension(vecLen*NSPECIES) :: massFrac
  ! must correspond to dimensions of Eos_wrapped
  logical,optional, dimension(EOS_VARS+1:EOS_NUM),INTENT(in)::mask

  integer :: i, k
  integer :: vecBegin,vecEnd
  integer :: pres, temp, dens, gamc, eint, game
  integer :: abar, zbar
  integer :: entr, dst, dsd
  integer :: dpt, dpd, det, ded, dea, dez, pel, ne, eta, c_v, c_p
  real    :: abarInv, zbarFrac
  real    :: xAbar, xZbar

  real, parameter :: KtoMeV = 1./temp_mev_to_kelvin
  real, parameter :: MeVtoK = temp_mev_to_kelvin

  real :: xDens,xTemp,xYe,xEner,xPres,xEntr,xCs2,xdedt, xGamc
  real :: xdpderho,xdpdrhoe,xmunu, xGame
  integer :: xMode, err
  real :: preTemp, chit, chid, x7
  logical :: doNuc

  vecBegin = 1
  vecEnd = vecLen

  err = 0

  ! These integers are indexes into the lowest location in UNK that contain the appropriate variable
  pres = (EOS_PRES-1)*vecLen
  dens = (EOS_DENS-1)*vecLen
  temp = (EOS_TEMP-1)*vecLen
  eint = (EOS_EINT-1)*vecLen   
  gamc = (EOS_GAMC-1)*vecLen   
  abar = (EOS_ABAR-1)*vecLen   
  zbar = (EOS_ZBAR-1)*vecLen   
  entr = (EOS_ENTR-1)*vecLen

  !      Crank the EOS on the pipes filled above, then fill the FLASH arrays
  !      with the thermodynamic quantities returned by the EOS.

  do k = vecBegin, vecEnd
! Get appropriate Ye
#ifdef YE_MSCALAR
     xYe = eosData(zbar+k) / eosData(abar+k) ! These are then set in Eos_getData
     xAbar = eosData(abar+k)
     xZbar = eosData(zbar+k)
#else
#ifdef FLASH_MULTISPECIES
     !Calculate the inverse in a way that allows for zero mass fractions
     call Multispecies_getSumInv(A, abarInv,massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     xAbar = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z, zbarFrac, massFrac((k-1)*NSPECIES+1:k*NSPECIES))
     xZbar = xAbar * zbarFrac
     xYe = xZbar / xAbar
#else
     ! No multispecies defined, use default values (same as Gamma formulation)
     xAbar = eos_singleSpeciesA
     xZbar = eos_singleSpeciesZ
     xYe = xZbar / xAbar 
#endif
#endif

     xDens = eosData(dens+k)
     xTemp = max(eosData(temp+k),eos_smallt)
     xEner = eosData(eint+k)
     xEntr = eosData(entr+k)
     xPres = eosData(pres+k)
     preTemp = xTemp

     doNuc = .false.
     select case(mode)
     case(MODE_DENS_EI)
        xMode = 0
     case(MODE_DENS_TEMP)
        xMode = 1
     case(MODE_DENS_ENTR)
        xMode = 2
        doNuc = .true. ! better know what you are doing when you call in this mode!
     case(MODE_DENS_PRES)
        xMode=4
     case default
        call Driver_abortFlash('[Eos] Error: unsupported mode for Nuclear Eos')
     end select

     ! Let's adjust the energy zero-point
     xEner = xEner - e_zeroPoint

     xTemp = xTemp * KtoMeV ! convert to MeV
     call nuc_eos_short(xDens,xTemp,xYe,xEner,xPres,xEntr,xCs2,xdedt,&
          xdpderho,xdpdrhoe,xmunu,xMode,err,precision)
     
     xGamc = xCs2*xDens/xPres
     
     xGame = xPres / (xDens*xEner) + 1.

     xTemp = xTemp * MeVtoK ! convert back to K

     if (err /= 0) then
        if (mode == MODE_DENS_ENTR .OR. mode == MODE_DENS_PRES) then
           ! For now, not aborting if an error comes back when
           ! when using dens+entr/pres.  Not sure of the side effects,
           ! but this seems to be fine, empirically speaking.
           ! We shall return from here so as to not accidentally
           ! change any in-coming eosData values.
           return
        else
           call Driver_abortFlash("[EOS] problem with nuclear EOS")
        endif
     endif
        
     xEner = xEner + e_zeroPoint ! re-adjust energy

     eosData(temp+k) = max(xTemp, eos_smallt)
     eosData(pres+k) = xPres
     eosData(eint+k) = xEner
     eosData(gamc+k) = xGamc
     eosData(entr+k) = xEntr
     eosData(abar+k) = xAbar
     eosData(zbar+k) = xZbar

  enddo

  return
end subroutine eos_nuclear
