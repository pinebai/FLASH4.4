!!****if* source/physics/materialProperties/Viscosity/ViscosityMain/Spitzer/Viscosity
!!
!!  NAME    
!!   Viscosity
!!
!!  SYNOPSIS
!!   call Viscosity(real(in) :: xtemp,
!!                  real(in) :: xden,
!!                  real(in) :: massfrac(NSPECIES),
!!                  real(out) :: viscDynamic,
!!                  real(out) :: viscKinematic)
!!                  
!!
!!  DESCRIPTION
!!   given temperature and density, this routine returns the
!!   Thermal viscosity given by Spitzer (1962)
!!
!!  ARGUMENTS
!!   xtemp    -  temperature (in K)
!!   xden     -  density (in g/cm**3)
!!   massfrac -   mass fractions of the composition
!!   viscDynamic   - returns coefficient of dynamic viscosity
!!   viscKinematic - returns coefficient of kinematic viscosity
!!
!!  NOTES
!!   See: Spitzer L. 1962, In: `Physics of fully ionized gases', 
!!   (New York: Wiley Interscience)
!!
!!  MODIFICATION HISTORY
!!   Written by: S. ORLANDO October 2001
!!
!!***

subroutine Viscosity(xtemp,xden,massfrac,viscDynamic,viscKinematic)
  
  use Viscosity_data, ONLY:  visc_useViscosity, viscSuppressFactor, &
       viscTempLow, viscTempHigh

  implicit none

#include "Flash.h"


  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac(NSPECIES)
  real,INTENT(out)   :: viscDynamic, viscKinematic
  
  real, PARAMETER :: ck = 1.25e-16 
  real, PARAMETER :: cexp = 2.5
  

  if (visc_useViscosity) then
     if (xtemp < viscTempLow .or. xtemp > viscTempHigh) then
        viscDynamic = 0.0
        viscKinematic = 0.0
     else
        viscDynamic   = viscSuppressFactor*ck*xtemp**cexp
        viscKinematic = viscDynamic / xden
     endif        
  else
     viscDynamic = 0.0
     viscKinematic = 0.0
  end if

end subroutine Viscosity
