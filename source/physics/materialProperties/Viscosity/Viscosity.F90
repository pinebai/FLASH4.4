!!****f* source/physics/materialProperties/Viscosity/Viscosity
!!
!!  NAME    
!!   Viscosity
!!
!!  SYNOPSIS
!!   Viscosity(real, intent(IN)  :: xtemp, 
!!             real, intent(IN)  :: xden, 
!!             real, inteng(IN)  :: massfrac(NSPECIES),
!!             real, intent(OUT) :: viscDynamic,
!!             real, intent(OUT) :: viscKinematic)
!!
!!
!! DESCRIPTION
!!   A generic viscosity routine.  Returns kinematic viscosity in
!!   viscKinematic, and the corresponding dynamic viscosity in viscDynamic.
!!
!! ARGUMENTS
!!
!!  INPUTS
!!   xtemp        :  temperature (in K)
!!   xden         :  density (in g/cm**3)
!!   massfrac     :  mass fractions of the composition
!!
!!  OUTPUTS
!!   viscDynamic  :  dynamic viscosity
!!   viscKinematic:  kinematic viscosity
!!
!!***

subroutine Viscosity(xtemp,xden,massfrac,viscDynamic,viscKinematic)
!! True Stub
  implicit none
#include "Flash.h"

  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac(NSPECIES)
  real,INTENT(out)   :: viscDynamic, viscKinematic

  !dummy value assigned in stub
  viscDynamic = 0.
  viscKinematic = 0.

  return 
end subroutine Viscosity
