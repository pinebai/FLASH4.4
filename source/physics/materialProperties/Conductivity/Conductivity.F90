!!****f* source/physics/materialProperties/Conductivity/Conductivity
!!
!! NAME
!!  Conductivity
!!
!! SYNOPSIS
!!  call Conductivity(real(in)    :: xtemp,
!!                    real(in)    :: xden,
!!                    real(in)    :: massfrac(NSPECIES),
!!                    real(out)   :: isochoricCond,
!!                    real(out)   :: diff_coeff,
!!                    integer(in) :: component)
!!
!! DESCRIPTION
!!
!!  Returns thermal conductivity and diffusivity coefficients.
!!
!!  The stub implementation for Conductivity returns isochoricCond = 0
!!  and diff_coeff = 0.
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond  :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)
!! true stub

#include "Flash.h"
  
  implicit none

  real, intent(IN)   :: xtemp
  real, intent(IN)   :: xden
  real, intent(IN)   :: massfrac(NSPECIES)
  real, intent(OUT)  :: diff_coeff
  real, intent(OUT)  :: isochoricCond
  integer,intent(IN) :: component


  !dummy values assigned in stub
  diff_coeff = 0.
  isochoricCond = 0.


  return 
end subroutine Conductivity

subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity

  implicit none

  real,              intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  real :: isochoricCondLoc, diffCoeffLoc
  integer :: componentLoc

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0
  
#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (present(component)) then
     componentLoc = component
  else
     componentLoc = 0
  end if
  call Conductivity(solnVec(TEMP_VAR),solnVec(DENS_VAR), &
       solnVec(SPECIES_BEGIN:SPECIES_END), &
       isochoricCondLoc, diffCoeffLoc, componentLoc)

#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc

end subroutine Conductivity_fullState

