!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity_fullState
!!
!! NAME
!!  Conductivity_fullState
!!
!! SYNOPSIS
!!  call Conductivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                     OPTIONAL,real(out)   :: isochoricCond,
!!                     OPTIONAL,real(out)   :: diffCoeff,
!!                     OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron conductivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!!  Returns thermal conductivity and/or diffusivity coefficients.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricCond  :   isochoric conductivity
!!   diffCoeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Flash.h"
#include "constants.h"  
#include "Eos.h"


subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar, Eos_getTempData

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  real, pointer :: massfrac(:)

  real :: isochoricCondLoc, diffCoeffLoc
  integer :: componentLoc, tempToUse

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen

  real :: xtemp, xden
  real :: nele, ll
  real :: abar, zbar
  real :: mion
  logical :: useTion, useTele
  integer :: tionVar, teleVar

  real, parameter :: cexp = 2.5

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0

#ifndef FLASH_3T
  call Driver_abortFlash("[Conductivity_fullState] SpitzerHighZ conductivity only works in 3T")
#endif

#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (present(component)) then
     componentLoc = component
  else
     componentLoc = 0
  end if

  tempToUse = TEMP_VAR

  select case (componentLoc)
  case(1)
#ifdef TION_VAR
     tempToUse = TION_VAR
#endif
  case(2)
#ifdef TELE_VAR
     tempToUse = TELE_VAR
#endif
  case(3)
#ifdef TRAD_VAR
     call Driver_abortFlash("[Conductivity_fullState] Spitzer conductivity does not work with radiation")
#endif
  end select

#ifdef TION_VAR
  useTion = (tempToUse == TION_VAR)
  tionVar = TION_VAR
#else
  useTion = .FALSE.
  tionVar = TEMP_VAR
#endif
#ifdef TELE_VAR
  useTele = (tempToUse == TELE_VAR)
  teleVar = TELE_VAR
#else
  useTele = .FALSE.
  teleVar = TEMP_VAR
#endif

    if (cond_useConductivity) then
       call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)
       xden = solnVec(DENS_VAR)
       nele = zbar * xden * cond_navo / abar

       if(useTele .or. tempToUse == TEMP_VAR) then

          call cond_logLambda(solnVec(teleVar), nele, zbar, ll)

          isochoricCondLoc = (8.0/PI)**1.5*cond_boltz**3.5 / (sqrt(cond_mele)*cond_qele**4) * &
               solnVec(teleVar)**cexp / (ll * (zbar + 3.3))

          if (present(diffCoeff)) then

             eos_arr(EOS_DENS) = xden
             call Eos_getTempData(solnVec,eos_arr,MODE_DENS_TEMP_GATHER)

             mask = .false.
             mask(EOS_CVELE)  = .true.
             mask(EOS_DET) = .true.

             if (NSPECIES > 0) then
                massfrac => solnVec(SPECIES_BEGIN:SPECIES_END)
                call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
             else
                call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,mask=mask)
             end if

             diffCoeffLoc = isochoricCondLoc/(xden*eos_arr(EOS_CVELE))
          end if

       elseif(useTion) then
          mion = abar / cond_navo

          ! Compute the electron conductivity:
          call cond_logLambdaII(solnVec(tionVar), solnVec(teleVar), nele, zbar, ll)
          isochoricCondLoc = &
               0.164 * 20.0 * (8.0/PI)**1.5*cond_boltz**3.5 * solnVec(tionVar)**2.5 / &
               (sqrt(mion) * cond_qele**4 * zbar**4 * ll) 

          if (present(diffCoeff)) then
             eos_arr(EOS_DENS) = xden
             call Eos_getTempData(solnVec,eos_arr,MODE_DENS_TEMP_GATHER)

             mask = .false.
             mask(EOS_CVION)  = .true.
             mask(EOS_DET) = .true.

             if (NSPECIES > 0) then
                massfrac => solnVec(SPECIES_BEGIN:SPECIES_END)
                call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
             else
                call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,mask=mask)
             end if

             diffCoeffLoc = isochoricCondLoc/(xden*eos_arr(EOS_CVION))
          end if
       end if
    end if

#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc
  
end subroutine Conductivity_fullState
