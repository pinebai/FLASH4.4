!!****if* source/physics/materialProperties/Viscosity/ViscosityMain/Constant/Viscosity
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
!!   A generic viscosity routine.  Just returns a constant kinematic viscosity
!!   from run-time parameter 'diff_visc_nu' in viscKinematic, and the
!!   corresponding dynamic viscosity in viscDynamic, if runtime parameter
!!   visc_whichCoefficientIsConst is 2.
!!   If visc_whichCoefficientIsConst is 1, use 'diff_visc_mu' as the constant
!!   dynamic viscosity instead, and compute a kinematic viscosity accordingly.
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

  use Viscosity_data, ONLY: visc_useViscosity, visc_diffNu,visc_diffMu, visc_whichCoefficientIsConst
  
  implicit none

#include "Flash.h"


  real,INTENT(out)   :: viscDynamic, viscKinematic
  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac(NSPECIES)

  if (visc_useViscosity) then

     ! Return constant viscosity
     if (visc_whichCoefficientIsConst==2) then
        viscKinematic = visc_diffNu
        viscDynamic   = visc_diffNu * xden
     else
        viscKinematic = visc_diffMu / xden
        viscDynamic   = visc_diffMu
     end if

  else
     ! Use of Viscosity is turned off
     viscDynamic = 0.0
     viscKinematic = 0.0
  end if

end subroutine Viscosity
