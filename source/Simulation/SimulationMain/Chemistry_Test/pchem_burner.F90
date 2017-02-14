!!***if* source/physics/sourceTerms/PrimordialChemistry/GlobChem/pchem_burner
!!
!! NAME
!!
!!  pchem_burner
!!
!! SYNOPSIS
!!
!!  pchem_burner(
!! real(in)   :: tstep,
!! real(in)   :: temp,
!!       real(in)   :: density,
!!       real(in)   :: xIn(:),
!! real(out)  :: xOut(:),
!! real(out)  :: sdotRate)
!!
!! DESCRIPTION
!!
!!  Routine pchem_burner drives the chemistry network
!!     given time step tstep, temperature temp, density denstiy, and
!!     composition xIn, this routine return the "chemed" composition xOut
!!     and the energy generation rate sdotRate.
!!
!! ARGUMENTS
!!
!! tstep:time step
!! temp:temperature
!! density: density
!! xIn:composition in
!! xOut: composition out
!! sdotRate:energy generation rate
!!
!! PARAMETERS
!!
!!  pchem_algebraintegerspecifies choice of pchem_algebra ma28=1 or gift=2
!!  odeStepper  integer specifies integration method bader-deuflhard=1 or rosenbrock=2
!!  useTableinteger NOT USED FOR THE CHEMISTRY I HAVE SET UP
!!
!!  NOTES
!!
!!  Within the network setup process
!!  routine pchem_network sets up the odes/derivatives
!!  routine pchem_networkDenseJackob sets up the dense jacobian
!!  routine pchem_networkRates sets up the raw reaction rates
!!  routine pchem_networkScreen applies screening corrections to the raw rates (none used)
!!  
!!  Within the general network integration process, found in directory ChemIntegrate:
!!    --> This is a direct copy of BurnIntegrate
!!  routine pchem_netIntegrate drives the integration of the odes
!!  this subroutine uses 8 routines called "steper" to drive the inegration of
!!  the reaction networks with
!!       2 choices of the time stepper and
!! 2 choices for the linear algebra (well not for us)
!!  routine pchem_baderMa28 drives a Bader-Deuflhard step with ma28 algebra
!!  routine pchem_baderStepMa28 is a Bader-Deuflhard stepper with ma28 algebra
!!  routine pchem_baderGift drives a Bader-Deuflhard step with gift algebra
!!  routine pchem_baderStepGift is a Bader_Deuflhard stepper with gift algebra
!!  routine pchem_rosenMa28 drives a Rosenbrock step with ma28 algebra
!!  routine pchem_rosenGift drives a Rosenbrock step with gift algebra
!!  FOR CHEMISTRY WE ARE LIMITED TO GIFT ALGEBRA (and trying Bader stepper)
!!  In addition,
!!  routine pchem_pzExtr does extrapolations for any of the Bader_Deuflhard pchem_bader* routines
!!  
!!  In this GlobChem directory, there are additional routines
!!  routine pchem_azbar computes composition variables from the mass fractions; they are stored
!!                in PrimordialChemistry_dataEOS
!!
!!***

subroutine pchem_burner(tstep,temp,density,xIn,xOut,sdotRate,ei,counts,jcounts,mfrac,tback)

  use PrimordialChemistry_dataEOS
  use PrimordialChemistry_data
  use Simulation_data
  use pchem_integrateInterface
  use pchem_networkInterface
  use pchem_interface
  use Eos_interface, ONLY: Eos
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
  !!external pchem_networkDenseJakob

  !arguments

  real, intent(IN)  :: tstep,temp,density
  real, intent(OUT)  :: sdotRate
  real, intent(IN), dimension(NSPECIES)   :: xIn
  real, intent(OUT), dimension(NSPECIES)  :: xOut
  real, intent(INOUT)  :: ei
  real, intent(OUT)  :: tback
  real, intent(INOUT)                     :: jcounts, counts, mfrac

  real   :: tfun

  !..local variables 
  integer:: i,k,nok,nbad,kount, NS
  integer, parameter :: tdim=10, iprint=0, nostore=0
  integer:: issmall

  real:: stptry,stpmin,ys2(NSPECIES),ttime(tdim),elem(NSPECIES,tdim)
  real, parameter       :: avo = 6.0221367e23, ev2erg = 1.602e-12, &
           conv = ev2erg*1.0e6*avo, tol = 1.0E-5, &
   beg = 0.0, odescal = 1.0E-6
  real:: temp_save, ei_save

  real:: ttt, yy(NSPECIES), ddyddx(NSPECIES)  !!DUMMY variables for use in error control 
  real:: dtprime, dtmax, ttmin, ttemp, smallx
  real                  :: eosData(EOS_NUM), xmsave(NSPECIES) , ei_out, ysav(NSPECIES), temp_back
  real:: xeos(SPECIES_BEGIN:SPECIES_END) !!For EOS calls
  real:: nden
  real:: eidot, tauei, eis
  real:: coolden


  !  print * , 'Here in PrimordialChemistry'
  !   print * , 'mfrac: ', mfrac , 'doCool: ', pchem_doCool 
  
  !..set the material and network variables
  tback = 0.0e0 !!send temperature back
  ctemp = temp
  den = density
  coolden = density
  denrate = density
  smallx = 1.0e-20
  !    ctemp= sim_c_temp
  !    den= sim_c_den
  
  !  print *,'ChemBurn temp: ', ctemp
  
  eosData(EOS_TEMP) = ctemp
  eosData(EOS_DENS) = den
  eosData(EOS_EINT) = ei
  
  !  print *, 'eosData(EOS_TEMP): ', temp , '    eosData(EOS_DENS): ', den,'    eosData(EOS_EINT): ', ei
  
  do i=1,NSPECIES
     xmass(i) = xIn(i)
     !    print *, 'xmass(',i,'): ' , xIn(i)
  enddo
  
  call pchem_azbar() !!get ymass
  
  do i=1,NSPECIES
     ys2(i) = ymass(i)
     if(ys2(i) .lt. 0.0) then
        ys2(i) = smallx
     endif
  enddo
  ys2(iELEC) = ys2(iHP)+ys2(iDP)+ys2(iHDP)+ys2(iH2P)+ys2(iHEP)+2.0*ys2(iHEPP)-ys2(iHM)-ys2(iDM)
  
  
  
  !..set the time step variables for a single point chem burn
  stptry = tstep
  stpmin = stptry*1.0e-20    !  stptry*1.0e-20
  
  ttt = 0.0
  counts = 0.0
  jcounts = 0.0
  
  if(ctemp .gt. 1.0E5) then
     counts = -1
     do i=1,NSPECIES
        ysav(i) = ys2(i)
     enddo
     ys2(iH) = smallx
     ys2(iHP) = sim_xH/aion(iH)
     ys2(iHM) = smallx
     ys2(iD) = smallx
     ys2(iDP) = sim_xD/aion(iD)
     ys2(iDM) = smallx
     ys2(iHE) = smallx
     ys2(iHEP) = smallx
     ys2(iHEPP) = sim_xHE/aion(iHEPP)
     ys2(iHD)  = smallx
     ys2(iHDP) = smallx
     ys2(iH2) = smallx
     ys2(iH2P) = smallx
     ys2(iELEC) = ys2(iHP) + ys2(iDP) + 2.0*ys2(iHEPP)
     sdotRate = 0.0e0
     do k=1,NSPECIES
        sdotRate = sdotRate + (ys2(k)-ysav(k))*bion(k)
     enddo
     xeos(H_SPEC) = ys2(iH)*aion(iH)
     xeos(HP_SPEC) = ys2(iHP)*aion(iHP)
     xeos(HM_SPEC) = ys2(iHM)*aion(iHM)
     xeos(D_SPEC) = ys2(iD)*aion(iD)
     xeos(DM_SPEC) = ys2(iDM)*aion(iDM)
     xeos(DP_SPEC) = ys2(iDP)*aion(iDP)
     xeos(HE_SPEC) = ys2(iHE)*aion(iHE)
     xeos(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
     xeos(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
     xeos(H2_SPEC) = ys2(iH2)*aion(iH2)
     xeos(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
     xeos(HD_SPEC) = ys2(iHD)*aion(iHD)
     xeos(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
     xeos(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
     
     ei = ei + sdotRate*conv
     eosData(EOS_EINT) = ei
     eosData(EOS_TEMP) = ctemp
     eosData(EOS_DENS) = den
     call Eos(MODE_DENS_TEMP,1,eosData,xeos)
     ctemp = eosData(EOS_TEMP)
     ei = eosData(EOS_EINT)
     den = eosData(EOS_DENS)  
     
     if(ctemp .lt. 50.0) then
        ctemp = 51.0
        eosData(EOS_TEMP) = ctemp
        eosData(EOS_DENS) = den
        eosData(EOS_EINT) = ei
        call Eos(MODE_DENS_TEMP,1,eosData,xeos)
        ei = eosData(EOS_EINT)
        den = eosData(EOS_DENS)
        ctemp = eosData(EOS_TEMP)
     endif
     
     
     coolden = den
     if(pchem_doCool .eq. 1) then
        call pchem_coolFunction(ctemp,coolden,ys2,ei,temp_back,ei_out,tstep,mfrac) !! GIVES NEW TEMP AND EI
        ei = ei_out
        ctemp = temp_back
     endif
     
     if(ctemp .lt. 50.0) then
        ctemp = 51.0
     endif
     
     eosData(EOS_EINT) = ei
     eosData(EOS_DENS) = den
     eosData(EOS_TEMP) = ctemp
     call Eos(MODE_DENS_TEMP,1,eosData,xeos)
     
     ei = eosData(EOS_EINT)
     ctemp = eosData(EOS_TEMP)
     den = eosData(EOS_DENS)
     
  else if(ctemp .lt. 50.0) then
     counts = -2
     ctemp = 51.0
     do k=1,NSPECIES
        if(ys2(k) .lt. 0.0) then
           ys2(k) = smallx
        endif
     enddo
     
     xeos(H_SPEC) = ys2(iH)*aion(iH)
     xeos(HP_SPEC) = ys2(iHP)*aion(iHP)
     xeos(HM_SPEC) = ys2(iHM)*aion(iHM)
     xeos(D_SPEC) = ys2(iD)*aion(iD)
     xeos(DM_SPEC) = ys2(iDM)*aion(iDM)
     xeos(DP_SPEC) = ys2(iDP)*aion(iDP)
     xeos(HE_SPEC) = ys2(iHE)*aion(iHE)
     xeos(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
     xeos(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
     xeos(H2_SPEC) = ys2(iH2)*aion(iH2)
     xeos(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
     xeos(HD_SPEC) = ys2(iHD)*aion(iHD)
     xeos(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
     xeos(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
     
     eosData(EOS_DENS) = den
     eosData(EOS_EINT) = ei
     eosData(EOS_TEMP) = ctemp
     !   print *, 'EOS call: COLD'
     call Eos(MODE_DENS_TEMP,1,eosData,xeos)
     ctemp = eosData(EOS_TEMP)
     ei = eosData(EOS_EINT)
     den = eosData(EOS_DENS)  
     
  else
     counts = 1.0
     
     if( ctemp .lt. 1.0E5 .and. ctemp .gt. 20000.0) then
        do i=1,NSPECIES    !!SAVE SPECIES
           ysav(i) = ys2(i)
        enddo
        
        if(ys2(iHP) .lt. 0.05*ys2(iH)) then
           ys2(iHP) = ys2(iHP) + 0.05*ys2(iH)
           ys2(iH)  = 0.95*ys2(iH)
        endif
        
        if((ys2(iHEPP)+ys2(iHEP)) .lt. 0.05*ys2(iHE)) then
           ys2(iHEP) = ys2(iHEP) +  0.045*ys2(iHE)
           ys2(iHEPP) = ys2(iHEPP)+ 0.005*ys2(iHE)
           ys2(iHE)   = 0.95*ys2(iHE)
        endif
        
        ys2(iELEC) = ys2(iHP) + ys2(iDP) + ys2(iHEP) + 2.0*ys2(iHEPP) + ys2(iHDP) + ys2(iH2P) - ys2(iHM) - ys2(iDM)
        
        sdotRate = 0.0e0
        do k=1,NSPECIES
           sdotRate = sdotRate + (ys2(k)-ysav(k))*bion(k)
        enddo
        
        xeos(H_SPEC) = ys2(iH)*aion(iH)
        xeos(HP_SPEC) = ys2(iHP)*aion(iHP)
        xeos(HM_SPEC) = ys2(iHM)*aion(iHM)
        xeos(D_SPEC) = ys2(iD)*aion(iD)
        xeos(DM_SPEC) = ys2(iDM)*aion(iDM)
        xeos(DP_SPEC) = ys2(iDP)*aion(iDP)
        xeos(HE_SPEC) = ys2(iHE)*aion(iHE)
        xeos(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
        xeos(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
        xeos(H2_SPEC) = ys2(iH2)*aion(iH2)
        xeos(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
        xeos(HD_SPEC) = ys2(iHD)*aion(iHD)
        xeos(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
        xeos(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
        
        ei = ei + sdotRate*conv
        eosData(EOS_EINT) = ei
        eosData(EOS_DENS) = den
        eosData(EOS_TEMP) = ctemp
        call Eos(MODE_DENS_TEMP,1,eosData,xeos)
        ctemp = eosData(EOS_TEMP)
        den = eosData(EOS_DENS)
        ei = eosData(EOS_EINT)
        
        if(ctemp .lt. 50.0) then
           ctemp = 51.0
           eosData(EOS_TEMP) = ctemp
           eosData(EOS_DENS) = den
           eosData(EOS_EINT) = ei
           call Eos(MODE_DENS_TEMP,1,eosData,xeos)
           ei = eosData(EOS_EINT)
           ctemp = eosData(EOS_TEMP)
           den = eosData(EOS_DENS)
        endif
        
     endif
     !!What the above does: If the temperature is between 2.0E4 and 1.0E5, then we check to see how much of the hydrogen
     !!and helium are ionized. If they are not ionized very much, then we ionize a little bit (5%)  and let PrimordialChemistry figure the 
     !!rest out. Hopefully, this will help speed up the network a bit
     
     
     call pchem_networkRates()          !!get the reaction rates from our equations
     call pchem_networkScreen(ymass)    !!Ok, just renames ratraw-->ratdums since there is no screening
     call pchem_network(ttt,ys2,ddyddx) !! This grabs my odes and puts them in a dummy array for me and it works!
     
     dtprime = tstep  !!Setup the initial time step
     ttemp = tstep
     tfun  = tstep
     issmall = -1 
     !!Need to find the inital timestep
     do i=1,NSPECIES-1
        ttmin = abs(sim_pchem_time*(ys2(i)+0.1*ys2(iHP))/ddyddx(i))
        if( ttmin .lt. ttemp) then 
           if(ys2(i) .gt. 1.0E-15) then 
              !!Just said if abundance is less than smallx, I don't care what your time step is
              ttemp = ttmin
              issmall = i
           endif
        endif
     enddo
     
     eis = 0.0e0
     do i=1,NSPECIES-1
        eis = eis + ddyddx(i)*bion(i)
     enddo
     tauei = 0.1*abs(ei/(eis*conv))   
     !      print *, 'FIRST TAUEI: ', tauei
     if( tauei .lt. ttemp) then
        ttemp = tauei
     endif
     
     do while(dtprime .gt. 0.0)    
        if(dtprime .le. ttemp) then  !NO SUB
           !  print *, 'NO SUB'
           !      tstep = dtprime
           !      stpmin = 1.0e-20*tstep
           stptry = dtprime
           stpmin = dtprime*1.0e-20 
           !..bader_deuflhard integration
           !..LEAVING OUT ANYTHING TO DO WITH MA28 ALGEBRA PACKAGE
           counts = counts
           
           do k=1,NSPECIES
              ysav(k) = ys2(k)
              !     print *, 'ysav: ', ysav(k), 'k: ', k
           enddo
           
           ! Choose beteween the two steppers and then you integrate foward in time by the hydro timestep
           if (pchem_odeStepper .eq. 1) then
              if (pchem_algebra .eq. 1) then
                 call Driver_abortFlash('ERROR in pchem_burner: using wrong Algebra package. Use GIFT')
              else if (pchem_algebra .eq. 2) then
                 call pchem_netIntegrate(beg,stptry,stpmin,dtprime,ys2,tol,beg,nostore,ttime,elem, &
                      tdim,NSPECIES,tdim,NSPECIES,nok,nbad,kount,odescal,iprint, &
                      pchem_network,pchem_networkDenseJakob,pchem_networkSparsePointers,pchem_baderGift,jcounts)
                 do k=1,NSPECIES-1
                    if(ys2(k) .lt. 0.0) then
                       ys2(k) = smallx
                    endif
                 enddo
                 ys2(iELEC) = ys2(iHP)+ys2(iDP)+ys2(iHDP)+ys2(iH2P)+ys2(iHEP)+2.0*ys2(iHEPP)-ys2(iHM)-ys2(iDM)
              end if
           else if (pchem_odeStepper .eq. 2 ) then
              if (pchem_algebra .eq. 1 ) then
                 call Driver_abortFlash('Error in pchem_burner2: using wrong Algebra package. Use GIFT')
                 
              else if (pchem_algebra .eq. 2) then
                 call pchem_netIntegrate(beg,stptry,stpmin,dtprime,ys2,tol,beg,nostore,ttime, &
                      elem,tdim,NSPECIES,tdim,NSPECIES,nok,nbad,kount,odescal, &
                      iprint,pchem_network,pchem_networkDenseJakob,pchem_networkSparsePointers,pchem_rosenGift,jcounts)
                 do k=1,NSPECIES 
                    if(ys2(k) .lt. 0.0) then
                       ys2(k) = smallx/aion(k)
                    endif
                 enddo
                 ys2(iELEC) = ys2(iHP)+ys2(iDP)+ys2(iHDP)+ys2(iH2P)+ys2(iHEP)+2.0*ys2(iHEPP)-ys2(iHM)-ys2(iDM)
              end if
              
           end if
           
           
           !!Time for new stuff. Going to call cooling here.
           
           
           sdotRate = 0.0e0
           do k=1,NSPECIES
              sdotRate = sdotRate + (ys2(k)-ysav(k))*bion(k)
              !     print *, 'SPECIES: ', k, 'BION: ', bion(k)
           enddo
           xeos(H_SPEC) = ys2(iH)*aion(iH)
           xeos(HP_SPEC) = ys2(iHP)*aion(iHP)
           xeos(HM_SPEC) = ys2(iHM)*aion(iHM)
           xeos(D_SPEC) = ys2(iD)*aion(iD)
           xeos(DM_SPEC) = ys2(iDM)*aion(iDM)
           xeos(DP_SPEC) = ys2(iDP)*aion(iDP)
           xeos(HE_SPEC) = ys2(iHE)*aion(iHE)
           xeos(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
           xeos(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
           xeos(H2_SPEC) = ys2(iH2)*aion(iH2)
           xeos(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
           xeos(HD_SPEC) = ys2(iHD)*aion(iHD)
           xeos(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
           xeos(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
           !!   print *, 'EI: ', ei, 'SDOT', sdotRate*conv, 'RATIO: ', ei/(sdotRate*conv)
           ei = ei + sdotRate*conv
           
           eosData(EOS_DENS) = den
           eosData(EOS_EINT) = ei
           eosData(EOS_TEMP) = ctemp
           
           call Eos(MODE_DENS_TEMP,1,eosData,xeos) 
           
           ei = eosData(EOS_EINT)
           ctemp = eosData(EOS_TEMP)
           den = eosData(EOS_DENS)
           
           if(ctemp .lt. 50.0) then
              ctemp = 51.0
              eosData(EOS_TEMP) = ctemp
              eosData(EOS_EINT) = ei
              eosData(EOS_DENS) = den
              call Eos(MODE_DENS_TEMP,1,eosData,xeos)
              ei = eosData(EOS_EINT)
              ctemp = eosData(EOS_TEMP)
              den = eosData(EOS_DENS)
           endif
           
           coolden = den
           
           if(pchem_doCool .eq. 1) then
              call pchem_coolFunction(ctemp,coolden,ys2,ei,temp_back,ei_out, dtprime, mfrac) !! GIVES NEW TEMP AND EI
              ei = ei_out
              ctemp = temp_back
              !    print *, 'CTEMP 1: ', ctemp
           endif
           
           if(ctemp .lt. 50.0) then
              ctemp = 51.0
              eosData(EOS_TEMP) = ctemp
              eosData(EOS_EINT) = ei
              eosData(EOS_DENS) = den
              call Eos(MODE_DENS_TEMP,1,eosData,xeos)
              ei = eosData(EOS_EINT)
              ctemp = eosData(EOS_TEMP)
              den = eosData(EOS_DENS)
           endif
           
           dtprime = 0.0
        else
           counts = counts + 1.0
           
           !  print *, 'CYCLE'
           
           ! In this case dtprime gt ttemp, so we have to subcycle
           do k=1,NSPECIES
              ysav(k) = ys2(k)
           enddo
           
           ! Inside of the ODE steppers, the full timestep is now ttemp
           stpmin = ttemp*1.0e-20
           stptry = ttemp
           
           if (pchem_odeStepper .eq. 1) then
              if (pchem_algebra .eq. 1) then
                 call Driver_abortFlash('ERROR in pchem_burner: using wrong Algebra package. Use GIFT')
              else if (pchem_algebra .eq. 2) then
                 call pchem_netIntegrate(beg,stptry,stpmin,ttemp,ys2,tol,beg,nostore,ttime,elem, &
                      tdim,NSPECIES,tdim,NSPECIES,nok,nbad,kount,odescal,iprint, &
                      pchem_network,pchem_networkDenseJakob,pchem_networkSparsePointers,pchem_baderGift,jcounts)
                 do k=1,NSPECIES
                    if(ys2(k) .lt. 0.0) then
                       ys2(k) = smallx
                    endif
                 enddo
                 ys2(iELEC) = ys2(iHP)+ys2(iDP)+ys2(iHDP)+ys2(iH2P)+ys2(iHEP)+2.0*ys2(iHEPP)-ys2(iHM)-ys2(iDM)
              end if
           else if (pchem_odeStepper .eq. 2 ) then
              if (pchem_algebra .eq. 1 ) then
                 call Driver_abortFlash('Error in pchem_burner2: using wrong Algebra package. Use GIFT')
              else if (pchem_algebra .eq. 2) then
                 call pchem_netIntegrate(beg,stptry,stpmin,ttemp,ys2,tol,beg,nostore,ttime, &
                      elem,tdim,NSPECIES,tdim,NSPECIES,nok,nbad,kount,odescal, &
                      iprint,pchem_network,pchem_networkDenseJakob,pchem_networkSparsePointers,pchem_rosenGift,jcounts)
                 do k=1,NSPECIES
                    if(ys2(k) .lt. 0.0) then
                       ys2(k) = smallx
                    endif
                 enddo
                 ys2(iELEC) = ys2(iHP)+ys2(iDP)+ys2(iHDP)+ys2(iH2P)+ys2(iHEP)+2.0*ys2(iHEPP)-ys2(iHM)-ys2(iDM)
              end if
           end if
           
           ! At this point chemistry has moved foward by a time ttemp
           ! so we subtract this from dtprime and we implement cooling and sdot
           dtprime = dtprime - ttemp
           ttemp = dtprime   !!rezero ttemp
           
           ! add cooling here
           sdotRate = 0.0e0
           do k=1,NSPECIES
              sdotRate = sdotRate + (ys2(k)-ysav(k))*bion(k)
              !      print *, 'SPECIE: ', k, 'bion(k): ', bion(k)
           enddo
           xeos(H_SPEC) = ys2(iH)*aion(iH)
           xeos(HP_SPEC) = ys2(iHP)*aion(iHP)
           xeos(HM_SPEC) = ys2(iHM)*aion(iHM)
           xeos(D_SPEC) = ys2(iD)*aion(iD)
           xeos(DM_SPEC) = ys2(iDM)*aion(iDM)
           xeos(DP_SPEC) = ys2(iDP)*aion(iDP)
           xeos(HE_SPEC) = ys2(iHE)*aion(iHE)
           xeos(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
           xeos(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
           xeos(H2_SPEC) = ys2(iH2)*aion(iH2)
           xeos(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
           xeos(HD_SPEC) = ys2(iHD)*aion(iHD)
           xeos(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
           xeos(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
           
           !    print *, 'SDOT: ', sdotRate
           ei = ei + sdotRate*conv !! *(ttemp/tstep)--> should not need this is sdotrate is zeroed above
           eosData(EOS_EINT) = ei
           eosData(EOS_DENS) = den
           eosData(EOS_TEMP) = ctemp
           !    print *, 'CALL EOS: SUBCYCLE'
           call Eos(MODE_DENS_TEMP,1,eosData,xeos)
           ei = eosData(EOS_EINT)
           ctemp = eosData(EOS_TEMP)
           den = eosData(EOS_DENS)
           
           if(ctemp .lt. 50.0) then
              ctemp = 51.0
              eosData(EOS_DENS) = den
              eosData(EOS_EINT) = ei
              eosData(EOS_TEMP) = ctemp
              call Eos(MODE_DENS_TEMP,1,eosData,xeos)
              ei = eosData(EOS_EINT)
              ctemp = eosData(EOS_TEMP)
              den = eosData(EOS_DENS)
           endif
           
           coolden = den
           
           if(pchem_doCool .eq. 1) then
              call pchem_coolFunction(ctemp,coolden,ys2,ei,temp_back,ei_out,ttemp,mfrac) !! GIVES NEW TEMP AND E    
              ei = ei_out
              ctemp = temp_back
           endif
           
           ! This is just a check to make sure that temperature doesn't go below 1
           if(ctemp .lt. 50.0) then
              ctemp = 51.0
              eosData(EOS_TEMP) = ctemp
              eosData(EOS_DENS) = den
              eosData(EOS_EINT) = ei
              call Eos(MODE_DENS_TEMP,1,eosData,xeos)   
              ei = eosData(EOS_EINT)
              ctemp = eosData(EOS_TEMP)
              den = eosData(EOS_DENS)
           endif
           
           !    print *, 'issmall= ', issmall, 'y(issmall) =' , ys2(issmall), 'y(iHP)= ', ys2(iHP), 'dydt(issmall)=', ddyddx(issmall)
           
           call pchem_networkRates()          !!get the reaction rates from our equations
           call pchem_networkScreen(ymass)    !!Ok, just renames ratraw-->ratdums since there is no screening
           call pchem_network(ttt,ys2,ddyddx)
           
           do i=1,NSPECIES-1
              ttmin = abs(sim_pchem_time*(ys2(i)+0.1*ys2(iHP))/ddyddx(i))
              if( ttmin .lt. ttemp) then 
                 if(ys2(i) .gt. 1.0e-15) then !!smallx) then 
                    !!Just said if abundance is less than smallx, I don't care what your time step is
                    ttemp = ttmin
                    issmall = i
                 endif
              endif
           enddo
           
           eis = 0.0e0
           do i=1,NSPECIES-1
              eis = eis + ddyddx(i)*bion(i)
           enddo
           tauei = 0.1*abs(ei/(eis*conv))   
           !     print *, 'SECOND TAUEI: ', tauei
           if( tauei .lt. ttemp) then
              ttemp = tauei
           endif
           
        endif
     enddo
     ! This ends the implementation of chemistry and cooling, 
     ! below we just write it back and implement a floor on abundances
     
     
  endif !!TEMP ENDIF
  
  xoktot  = xoktot + real(nok)
  xbadtot = xbadtot + real(nbad)
  
  !Put in a floor for abundances
  do k = 1, NSPECIES
     if(ys2(k) .lt. 0.0) then
        ys2(k) = smallx
     endif
  enddo
  
  !   print *, 'COUNTS: ', counts
  
  !.. average energy generated over the time step
  sdotRate = 0.0e0
  do k=1,NSPECIES
     sdotRate = sdotRate + (ys2(k)-ymass(k))*bion(k)
  enddo
  
  sdotRate = sdotRate*conv/tstep
  
  !..Not going to worry right now about neutrino loss
  
  !..update composition
  do i=1,NSPECIES
     xOut(i) = ys2(i)*aion(i)  
  enddo
  
  tback = ctemp  
  !    xOut(H_SPEC) = ys2(iH)*aion(iH)
  !    xOut(HP_SPEC) = ys2(iHP)*aion(iHP)
  !    xOut(HM_SPEC) = ys2(iHM)*aion(iHM)
  !    xOut(D_SPEC) = ys2(iD)*aion(iD)
  !    xOut(DM_SPEC) = ys2(iDM)*aion(iDM)
  !    xOut(DP_SPEC) = ys2(iDP)*aion(iDP)
  !    xOut(HE_SPEC) = ys2(iHE)*aion(iHE)
  !    xOut(HEP_SPEC) = ys2(iHEP)*aion(iHEP)
  !    xOut(HEPP_SPEC) = ys2(iHEPP)*aion(iHEPP)
  !    xOut(H2_SPEC) = ys2(iH2)*aion(iH2)
  !    xOut(H2P_SPEC) = ys2(iH2P)*aion(iH2P)
  !    xOut(HD_SPEC) = ys2(iHD)*aion(iHD)
  !    xOut(HDP_SPEC) = ys2(iHDP)*aion(iHDP)
  !    xOut(ELEC_SPEC) = ys2(iELEC)*aion(iELEC)
  
  return
  
end subroutine pchem_burner

        
