!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant-diff/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  Conductivity(real, intent(IN)  :: xtemp,
!!               real, intent(IN)  :: xden,
!!               real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!               real, intent(OUT)  :: isochoricCond,
!!               real, intent(OUT)  :: diff_coeff,
!!               integer,intent(IN) :: component)
!!
!! DESCRIPTION
!!   
!!   Applies conductivity with constant coefficient of diffusivity
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond  :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)
  
  use Conductivity_data,ONLY :cond_diffConstant, cond_useConductivity
  use Eos_interface, ONLY : Eos
  implicit none
  
#include "constants.h"  
#include "Flash.h"
#include "Eos.h"
  
  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, isochoricCond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen

  if (cond_useConductivity) then
     diff_coeff = cond_diffConstant
     vecLen = 1
     mode = MODE_DENS_TEMP
     eos_arr(EOS_TEMP) = xtemp
     eos_arr(EOS_DENS) = xden

     mask = .false.
     mask(EOS_CV) = .true.
!!     mask(EOS_CP) = .true.
     mask(EOS_DET) = .true.
     call Eos(mode,vecLen,eos_arr,massfrac, mask)
     isochoricCond = (cond_diffConstant*eos_arr(EOS_CV)*xden)

  else
     isochoricCond = 0.0
     diff_coeff = 0.0
  end if

end subroutine Conductivity
