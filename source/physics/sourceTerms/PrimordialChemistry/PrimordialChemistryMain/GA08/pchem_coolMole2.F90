!!***if* source/physics/sourceTerms/PrimordialChemistry/PrimordialChemistryMain/GA08/pchem_coolMole2
!!
!! NAME
!!  pchem_coolMole2
!!
!!
!! SYNOPSIS
!!
!!  call pchem_coolMole2(temp, rho, yIn(ispecies),cool_out)
!! 
!!
!! 
!! DESCRIPTION
!!
!! this routine calulcates all the cooling rates and figures out the change in
!! energy, edited for coupling with chemistry 
!!
!!
!!
!!
!!
!!***

subroutine pchem_coolMole2(temp,rho,yIn,cool_out)

use PrimordialChemistry_data

implicit none

#include "constants.h"
#include "Flash.h"

         real   :: hnden, hpnden, hmnden, dnden, dpnden, dmnden, &
                   henden, hepnden,heppnden, hdnden, hdpnden, &
                   h2nden, h2pnden, enden
         real   :: cool_sum
         real   :: cool_term,cool_term_sum
         real   :: T3
         real   :: Na
         real   :: lt1, lt2, lt3, lt4, lt5
         real   :: temp,rho, cool_out
         real, dimension(NSPECIES) :: yIn

         Na = 6.0221415e23
         T3 = temp/1000.0

         cool_sum = 0.0
         cool_term = 0.0
         cool_term_sum = 0.0 
         lt1 = log10(T3)
         lt2 = log10(T3)**2
         lt3 = log10(T3)**3
         lt4 = log10(T3)**4
         lt5 = log10(T3)**5

        hnden    = yIn(iH)   * rho*Na
        hmnden   = yIn(iHM)  * rho*Na
        hpnden   = yIn(iHP)  * rho*Na
        dnden    = yIn(iD)   * rho*Na
        dmnden   = yIn(iDM)  * rho*Na
        dpnden   = yIn(iDP)  * rho*Na
        henden   = yIn(iHE)  * rho*Na
        hepnden  = yIn(iHEP) * rho*Na
        heppnden = yIn(iHEPP)* rho*Na
        hdnden   = yIn(iHD)  * rho*Na
        hdpnden  = yIn(iHDP) * rho*Na
        h2nden   = yIn(iH2)  * rho*Na
        h2pnden  = yIn(iH2P) * rho*Na
        enden    = yIn(iELEC)* rho*Na
        
! First: H2 and H

         if(temp .gt. 10.0 .and. temp .lt. 100.0) then
            cool_term = -16.818342 + 37.383713*lt1 + 58.145166*lt2 + 48.656103*lt3 + 20.159831*lt4 + 3.8479610*lt5
            cool_term_sum = cool_term_sum+10**(cool_term)
            cool_term = 10.0**(cool_term)*h2nden*hnden/rho
            cool_sum = cool_sum + cool_term
         else if(temp .gt. 100.0 .and. temp .lt. 1000.0) then
            cool_term = -24.311209 + 3.5692468*lt1 - 11.332860*lt2 -27.850082*lt3 -21.328264*lt4 -4.2519023*lt5
            cool_term_sum = cool_term_sum+10**(cool_term)
            cool_term = 10.0**(cool_term)*h2nden*hnden/rho
            cool_sum = cool_sum + cool_term
         else if(temp .gt. 1000.0 .and. temp .lt. 6000.0) then
            cool_term = -24.311209 + 4.6450521*lt1 -3.7209846*lt2 +5.9369081*lt3 -5.5108047*lt4 +1.5538288*lt5
            cool_term_sum = cool_term_sum+10**(cool_term)
            cool_term = 10.0**(cool_term)*h2nden*hnden/rho
            cool_sum = cool_sum+cool_term
         else
            cool_term = 0.0
            cool_sum = cool_term+cool_sum
            cool_term_sum = cool_term_sum
         endif

!Next: H2 and H2
         if(temp .gt. 100.0 .and. temp  .lt. 6000.0) then
            cool_term = -23.962112  + 2.09433740*lt1 -0.77151436*lt2 +0.43693353*lt3 -0.14913216*lt4 -0.033638326*lt5
            cool_term_sum = cool_term_sum+10**(cool_term)
            cool_term = 10**(cool_term)*h2nden*h2nden/rho
            cool_sum  = cool_sum + cool_term
         else
            cool_term = 0.0;
            cool_sum = cool_sum + cool_term
            cool_term_sum = cool_term_sum
         endif
         

!Next: H2 and He
       if(temp .gt. 10.0 .and. temp  .lt. 6000.0) then
            cool_term = -23.689237 + 2.1892372*lt1 -0.81520438*lt2 +0.29036281*lt3 -0.16596184*lt4 +0.19191375*lt5
            cool_term_sum = cool_term_sum+10**(cool_term)
            cool_term = 10.0**(cool_term)*h2nden*henden/rho
            cool_sum = cool_sum + cool_term
       else
            cool_term = 0.0;
            cool_sum = cool_sum + cool_term
            cool_term_sum = cool_term_sum
       endif

!Next: H2 and H+
       if(temp .gt. 10.0 .and. temp  .lt. 10000.0) then
         cool_term = -21.716699 +1.3865783*lt1 -0.37915285*lt2 +0.11453688*lt3 -0.23214154*lt4 +0.058538864*lt5
         cool_term_sum = cool_term_sum+10**(cool_term)
         cool_term = 10.0**(cool_term)*h2nden*hpnden/rho
         cool_sum = cool_sum + cool_term
       else
         cool_term = 0.0;
         cool_sum = cool_sum+cool_term
         cool_term_sum = cool_term_sum
       endif


!Next: H2 and e
       
         if( temp .gt. 10.0 .and. temp .lt. 200.0) then
           cool_term = -34.286155 -48.537163*lt1 -77.1251176*lt2 -51.352459*lt3 -15.169160*lt4 -0.98120322*lt5
           cool_term_sum = cool_term_sum+10**(cool_term)
           cool_term = 10.0**(cool_term)*h2nden*enden/rho
           cool_sum = cool_sum + cool_term
         else if(temp .gt. 200.0 .and. temp .lt. 10000.0) then
           cool_term = -22.190316 +1.5728955*lt1 -0.21335100*lt2 +0.96149759*lt3 -0.91023195*lt4 +0.13749749*lt5
           cool_term_sum = cool_term_sum+10**(cool_term)
           cool_term = 10.0**(cool_term)*h2nden*enden/rho
           cool_sum = cool_sum + cool_term
         else
           cool_sum = cool_sum
         endif

! Lipovka, Nunez, Reese Cooling rate for HD and H
        if(temp .gt. 30.0 .and. temp .lt. 20000.0) then
           cool_term = -42.45906 + 21.90083*log10(temp) -10.1954*log10(temp)**2 +2.19788*log10(temp)**3 -0.17286*log10(temp)**4
           cool_term_sum = cool_term_sum+10**(cool_term)
           cool_term = 10.0**(cool_term)*hdnden*hnden/rho
           cool_sum = cool_sum + cool_term
        else
           cool_sum = cool_sum
        endif

! Electron H2+
           if(temp .lt. 2000.0 .and. temp .gt. 1.0) then
                   cool_term = 1.1e-19*temp**(-0.34)*exp(-3025.0/temp)
                   cool_term_sum = cool_term_sum+(cool_term)
                   cool_term = (cool_term)*enden*h2pnden/rho
                    cool_sum = cool_sum + cool_term
           else
                cool_term = 3.35e-21*temp**(0.12)*exp(-3025.0/temp)
                cool_term_sum = cool_term_sum+(cool_term)
                cool_term = (cool_term)*enden*h2pnden/rho
                cool_sum = cool_sum + cool_term
           endif

! H and H2+
           
           if(temp .lt. 1000.0 .and. temp .gt. 1.0) then
                cool_term = 1.36e-22*exp(-3152.0/temp)
                cool_term_sum = cool_term_sum+(cool_term)
                cool_term = cool_term*hnden*h2pnden/rho
                cool_sum = cool_sum + cool_term
           else
                cool_term = 10.0**(-36.42 + 5.95*log10(temp) -0.526*log10(temp)**2);
                cool_term_sum = cool_term_sum+(cool_term)
                cool_term = cool_term*hnden*h2pnden/rho
                cool_sum = cool_sum + cool_term
           endif

         
         cool_out = cool_sum !!The real one to use!
!        cool_out = 0.0
!        cool_out = cool_term_sum  !!THIS IS ON JUST TO TEST COOLING FUNCTION NOT TOTAL COOLING

return
end subroutine pchem_coolMole2
