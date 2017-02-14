!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw-gray/Conductivity
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
!!   Thermal conductivity and diffusion coefficient given by a power law;
!!   this includes Spitzer (1962) by choosing runtime parameters appropriately.
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   cond       :   conductivity
!!   diff_coeff :   diffusion coefficient ( = cond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!  NOTES
!!   See: Spitzer L. 1962, In: `Physics of fully ionized gases', 
!!   (New York: Wiley Interscience)
!!
!!  MODIFICATION HISTORY
!!
!!   Initial: K. Weide March 2010
!!
!!***



subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff, component)

  use Conductivity_data, ONLY: cond_useConductivity, cond_meshMe, &
       cond_TemperatureExponent, cond_K0, cond_alpha, Raddiff_K0r, Raddiff_TemperatureExponent
  use Eos_interface, ONLY : Eos
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
 
  implicit none
  
#include "constants.h"  
#include "Flash.h"
#include "Eos.h"
  
  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen
  
  real :: ck, C, cv
  real :: cexp


  !Based on a switch we determine, which to compute.
  if (component == 3) then !! Rad Diffusivity

      ! Rad DiffCoeff = C / 3*K, C-Speed of light, K-Opacity, K=K0*T**n, 
      ! T would be the material temperature.
      ! K0,n can be specified in flash.par

      !call PhysicalConstants_get("speed of light",C)

      !cond = Raddiff_K0r*(xtemp**Raddiff_TemperatureExponent)

      !diff_coeff = (C/(3.0*cond))

      cexp = Raddiff_TemperatureExponent
      ck   = Raddiff_K0r
      cond       = ck*xtemp**cexp
      vecLen = 1
      mode = MODE_DENS_TEMP
      eos_arr(EOS_TEMP) = xtemp
      eos_arr(EOS_DENS) = xden

      mask = .false. 
      mask(EOS_CV) = .true.
      !! mask(EOS_CP) = .true.
      mask(EOS_DET) = .true.

      call Eos(mode,vecLen,eos_arr,massfrac,mask)
      diff_coeff = cond/(xden*eos_arr(EOS_CV))

     
  else if (component == 2) then      ! Electron Conductivity  (Power Law)

      cexp = cond_TemperatureExponent
      ck   = cond_K0
      cond       = ck*xtemp**cexp
      vecLen = 1
      mode = MODE_DENS_TEMP
      eos_arr(EOS_TEMP) = xtemp
      eos_arr(EOS_DENS) = xden

      mask = .false.
#ifdef TELE_VAR
      mask(EOS_DET) = .true.
      mask(EOS_CVELE) = .true.
      call Eos(mode,vecLen,eos_arr,massfrac,mask)      
      cv = eos_arr(EOS_CVELE)
#else
      mask(EOS_CV) = .true.
      mask(EOS_DET) = .true.
      call Eos(mode,vecLen,eos_arr,massfrac,mask)      
      cv = eos_arr(EOS_CV)
#endif

      diff_coeff = cond/(xden*cv)

   else

      if (cond_meshMe == MASTER_PE) &
           print*,'Conductivity: requested invalid component',component
      call Driver_abortFlash("Conductivity: requested invalid component")

 endif

end subroutine Conductivity
