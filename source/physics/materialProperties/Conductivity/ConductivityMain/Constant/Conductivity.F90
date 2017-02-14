!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant/Conductivity
!!
!! NAME
!!
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
!!   This routine applies constant conductivity
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond       :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!  NOTES
!!
!!
!!***

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)

  ! a generic conductivity routine.  Just return a constant isochoric conductivity.
  ! Note that isochoricCond is constant but the coefficient of of diffusivity, diff_coeff,
  ! is not (in general)!
  
  use Conductivity_data, ONLY: cond_constantIsochoric, cond_useConductivity
  use Eos_interface, ONLY : Eos
  implicit none

#include "constants.h"  
#include "Flash.h"
#include "Eos.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, isochoricCond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  real, dimension(EOS_NUM) :: eos_arr
  integer,intent(IN) :: component

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen
  
  if (cond_useConductivity) then
     vecLen = 1
     mode = MODE_DENS_TEMP
     eos_arr(EOS_TEMP) = xtemp
     eos_arr(EOS_DENS) = xden
     eos_arr(EOS_ABAR) = 1.0    !DEV: ???
     eos_arr(EOS_ZBAR) = 1.0    !DEV: ???

     
     mask = .false.
     mask(EOS_CV) = .true.
!!     mask(EOS_CP) = .true.
     mask(EOS_DET) = .true.

     call Eos(mode,vecLen,eos_arr,massfrac,mask)

     isochoricCond = cond_constantIsochoric
     diff_coeff = isochoricCond /(xden*eos_arr(EOS_CV))
  else
     isochoricCond = 0.0
     diff_coeff = 0.0
  end if
 

end subroutine Conductivity
