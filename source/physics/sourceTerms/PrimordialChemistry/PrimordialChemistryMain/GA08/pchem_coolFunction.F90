!!****if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_coolFunction
!!
!! NAME
!!  
!!  pchem_coolFunction
!!
!!
!! SYNOPSIS
!! 
!!  call pchem_coolFunction(temperture, density, y(NSPECIES), ei, temp_out, ei_out, times, mfrac)
!!  
!! DESCRIPTION
!!
!!  Apply cooling to a cell. There are 3 cooling terms, 1: High-temperature H-HE cooling (both case A and case B)
!!  2: Metal-line cooling via a METL_MSCALAR, and 3: Molecular line cooling
!!
!!  
!!
!!
!!  After we call radloss, call the eos to update the pressure 
!!  and temperature based on the radiative losses.
!!
!!
!!***

 subroutine pchem_coolFunction(tmp, rho, yin, ei ,temp_out, ei_out, times, mfrac)
      use PrimordialChemistry_data
      use Driver_interface, ONLY: Driver_getDt
      use Grid_interface, ONLY: Grid_getBlkIndexLimits, Grid_getBlkPtr, &
        Grid_releaseBlkPtr
      use Eos_interface, ONLY : Eos
      use Simulation_data, ONLY: sim_cool_time
      implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

      real ::  yin(NSPECIES)

      real :: radia,sdot,cool_out

!     These are counters
      integer :: i, j, k

!     Specie Stuff
       real                             :: xIn(SPECIES_BEGIN:SPECIES_END)
       integer :: n,specieMap

!     This is the time step
      real :: dt, t_save
      real :: dtprime
      real :: dtmax

      real :: tmp, rho, den_H, den_e, ei, eiold, ek, temp_out, ei_out, times, tmpsav, mfrac
      real :: indexfloat
      integer :: index

!     This is for low-temperature Metals
      EXTERNAL lowt_metals  !!The function
      real :: lowt_metals
      real :: ltm, lowmet          !!The return

! Mark zones that suffer cooling
      logical :: rad_zone
      real, dimension(EOS_NUM) :: eosData
      external radloss
!
!==============================================================================
!
! get the current timestep
     dt = times
     rad_zone = .FALSE.
     t_save = times
! sweep over all the zones



    
! Pass tmp and rho into function

                

                xIn(H_SPEC)    =        yin(iH)    * aion(iH)
                xIn(HM_SPEC)   =        yin(iHM)   * aion(iHM)
                xIn(HP_SPEC)   =        yin(iHP)   * aion(iHP)
                xIn(D_SPEC)    =        yin(iD)    * aion(iD)
                xIn(DM_SPEC)   =        yin(iDM)   * aion(iDM)
                xIn(DP_SPEC)   =        yin(iDP)   * aion(iDP)
                xIn(HE_SPEC)   =        yin(iHE)   * aion(iHE)
                xIn(HEP_SPEC)  =        yin(iHEP)  * aion(iHEP)
                xIn(HEPP_SPEC) =        yin(iHEPP) * aion(iHEPP)
                xIn(H2_SPEC)   =        yin(iH2)   * aion(iH2)
                xIn(H2P_SPEC)  =        yin(iH2P)  * aion(iH2P)
                xIn(HD_SPEC)   =        yin(iHD)   * aion(iHD)
                xIn(HDP_SPEC)  =        yin(iHDP)  * aion(iHDP)
                xIn(ELEC_SPEC) =        yin(iELEC) * aion(iELEC)

                eosData(EOS_TEMP) = tmp
                eosData(EOS_DENS) = rho


               if(tmp .lt. pchem_tradmin) then
                 if(tmp .lt. 50.0) then
                    tmp = 51.0
                    eosData(EOS_TEMP) = tmp
                    call Eos(MODE_DENS_TEMP,1,eosData,xIn)
                    ei = eosData(EOS_EINT)
                 endif  
                rad_zone = .TRUE.
               endif 

             !!  print *, 'RHO: ', rho, 'masfracH: ', pchem_massFracH, 'amu: ', amu

! derive the Hydrogen and the electron number density
              den_H = rho*pchem_massFracH/amu
              den_e = rho*(pchem_massFracH+(2.*(1.-pchem_massFracH)/4.))/amu 
             
              sdot = 0.0e0

! if the temperature is between the limits and
! the density is between the limits

             !!  print *, 'pchem_tradmin: ', pchem_tradmin, 'pchem_tradmax: ', pchem_tradmax
             !!  print *, 'den_h: ', den_h, 'amu: ', amu
             !!  print *, 'dmin: ', den_H*amu, ' dmax: ', den_H*amu
             !!  print *, 'pchem_dradmin: ', pchem_dradmin, 'pchem_dradmax: ', pchem_dradmax

               if ( (tmp >= pchem_tradmin .AND. tmp <= pchem_tradmax) .AND.         &
     &        (den_H*amu >= pchem_dradmin .AND. den_H*amu <= pchem_dradmax) ) then
              !! print *, 'HERE IN COOL_FUNC'
              !!   print *, 'pchem_noCool: ', pchem_noCool
                 rad_zone = .TRUE.

! radiative losses from an optically thin plasma
                  radia = 0.0
                  if(tmp>=pchem_tradmin) then
                    indexfloat = 353.*(log(tmp)/log(10.)-2.)/7.
                    index = indexfloat
                    if(index>=353) then 
                      index = 353
                    endif
                    if(index .lt. 1) then
                      index = 1
                    endif

                    if( tmp > 6000.0) then
                      if(pchem_ccCase .eq. 1) then
                        radia =(will_hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2)))
                      else if(pchem_ccCase .eq. 0) then
                         radia = hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2))
                      else
                         radia =(will_hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2)))
                      endif
                    else
                       if(tmp .gt. 100.0 .AND. tmp .lt. 6.0E3) then
                         radia = metal_cooling(index)*mfrac/(1.22E-2)
                       endif
                    endif
                  endif

                  sdot = -(den_H*den_e*radia)/rho

                 ! print *, 'COOLFUN: '
                 ! print *, 'radia: ', radia, 'sdot: ', sdot
! MY TURN
           
               if(tmp .gt. pchem_tradmin) then
!CALL MY FUNCTION HERE
                 call pchem_coolMole2(tmp,rho,yin,cool_out)
                 sdot = sdot - cool_out
               endif

! New stuff: This will loop over time to make sure we don't cool too fast and get wrong answers

                  dtprime = dt
                  dtmax = sim_cool_time*ei/abs(-sdot*pchem_noCool+1E-30)
                 do while(dtprime> 0.0) 
                     if(dtprime < dtmax) then
                        ei = ei + dtprime*sdot*pchem_noCool
                        eosData(EOS_EINT) = ei
                        eosData(EOS_TEMP) = tmp
                        eosData(EOS_DENS) = rho
                        call Eos(MODE_DENS_EI,1,eosData,xIn)
                        ei = eosData(EOS_EINT)
                        tmp = eosData(EOS_TEMP)
                        rho = eosData(EOS_DENS)
                        dtprime = 0.0
                     else
                        ei = (1.0-sim_cool_time)*ei
                        tmp = (1.0-sim_cool_time)*tmp
                        dtprime = dtprime-dtmax
                        den_H = rho*pchem_massFracH/amu
                        den_e = rho*(pchem_massFracH+(2.*(1.-pchem_massFracH)/4.))/amu
                        sdot = 0.0
                        radia = 0.0
                        if(tmp>=pchem_tradmin) then
                           indexfloat = 353.*(log(tmp)/log(10.)-2.)/7.
                           index = indexfloat
                          if(index>=353) then 
                              index = 353
                          endif
                          if(index .lt. 1) then
                            index = 1
                          endif

                           if( tmp > 6000.0) then
                             if(pchem_ccCase .eq. 1) then
                                radia = (will_hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2)))
                             else if(pchem_ccCase .eq. 0) then
                                radia = (hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2)))
                             else
                                radia = (will_hhe_cooling(index)+metal_cooling(index)*(mfrac/(1.22E-2)))
                             endif
                           else
                             if(tmp .gt. 100.0 .AND. tmp .lt. 6.0E3) then
                                radia = metal_cooling(index)*mfrac/(1.22E-2)
                             endif
                           endif
                         endif
                  sdot = -(den_H*den_e*radia)/rho
! MY TURN
                  if(tmp .gt. pchem_tradmin) then
!CALL MY FUNCTION HERE
                   call pchem_coolMole2(tmp,rho,yin,cool_out)
                   sdot = sdot - cool_out
                  endif

                   dtmax = sim_cool_time*ei/abs(-sdot*pchem_noCool+1.0E-30)
                  endif
                 enddo

               endif
      
       temp_out = tmp
       ei_out = ei

     return
end subroutine pchem_coolFunction
