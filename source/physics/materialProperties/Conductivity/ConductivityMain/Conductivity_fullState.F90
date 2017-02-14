!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Conductivity_fullState
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
!!  Returns thermal conductivity and/or diffusivity coefficients.
!!
!!  The stub implementation for Conductivity_fullState returns isochoricCond = 0
!!  and diffCoeff = 0.
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
!!
!!
!!***

#include "Flash.h"

subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity

  implicit none

  real,              intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  real :: isochoricCondLoc, diffCoeffLoc
  integer :: componentLoc
  integer :: tempToUse

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0
  tempToUse = TEMP_VAR
  
#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (present(component)) then
     componentLoc = component
  else
     componentLoc = 0
  end if

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
        tempToUse = TRAD_VAR
#endif
     end select

  call Conductivity(solnVec(tempToUse),solnVec(DENS_VAR), &
       solnVec(SPECIES_BEGIN:SPECIES_END), &
       isochoricCondLoc, diffCoeffLoc, componentLoc)

#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc

end subroutine Conductivity_fullState

