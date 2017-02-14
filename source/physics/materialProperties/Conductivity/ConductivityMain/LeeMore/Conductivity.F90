!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/LeeMore/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  call Conductivity(real, intent(IN)  :: xtemp,
!!                    real, intent(IN)  :: xden,
!!                    real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!                    real, intent(OUT)  :: cond,
!!                    real, intent(OUT)  :: diff_coeff,
!!                    integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Lee & More (Phys. Fluids 27, 1273 1984) electron thermal conductivity
!! for all materials. 
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
!!***

#include "Flash.h"

subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff, component)
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo, cond_hbar
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar
!  use Simulation_data, ONLY: sim_zbarTarg, sim_abarTarg
  use Multispecies_interface, ONLY: Multispecies_getSumInv, &
      Multispecies_getSumFrac

  implicit none
  
#include "constants.h"  
#include "Eos.h"
#include "Multispecies.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen

  real :: nele, ll
  real :: abar, zbar, znucbar, abarValue, abarInv, zbarFrac 
  real :: abarTarg, zbarTarg
  real :: tau, A_beta, nion, relax_time, R0
  real :: mu_div_kT, Fermi_onehalf, n3D, y
  real :: b, xi, Tmelt

  real, parameter :: a1 = 13.5
  real, parameter :: a2 = 0.976
  real, parameter :: a3 = 0.437
  real, parameter :: b2 = 0.51
  real, parameter :: b3 = 0.126

  if (cond_useConductivity) then
     call Eos_getAbarZbar(abar=abar,zbar=zbar,massFrac=massfrac)

     call Multispecies_getSumInv(A, abarInv ,massfrac)
     abarValue = 1.e0 / abarInv

     call Multispecies_getSumFrac(Z,zbarFrac,massfrac)
     znucbar = abarValue*zbarFrac

     abarTarg = abarValue
     zbarTarg = znucbar

     nele = zbar * xden * cond_navo / abar

     !Passing xtemp as both the electron and ion temperature
     call cond_logLambda(xtemp, nele, xtemp, zbar, ll)

     n3D = 2 * (cond_mele * cond_boltz * xtemp / (2 * PI * cond_hbar**2 ))**(3./2)
     Fermi_onehalf = nele/n3D

     call fermi_chem_potential(nele, xtemp, mu_div_kT)

     ! LM84, Eq. 27 including g(Z) = 1/(1+3.3/zbar) spitzer correction (c.f. Kruer, Atzeni books)
     tau = 3 /(2*sqrt(2.0)*PI) * cond_boltz**(1.5) * sqrt(cond_mele) / &
                (cond_qele**4) *xtemp**(1.5) / (ll * (zbar+3.3) * nele)  * &
                        (1 + exp(-mu_div_kT)) * Fermi_onehalf
 
     ! c.f. LM84 Appendix A
     if (mu_div_kT > 20) then
          y = mu_div_kT
     else if (mu_div_kT < -15) then
          y = exp(mu_div_kT)
     else 
          y = log(1.+exp(mu_div_kT));
     end if
     A_beta = (a1 + a2*y + a3*y**2 ) / (1 + b2*y + b3*y**2);  

     !Spitzer limit, LM84 Eq. 28a, for debugging
     !A_beta = 128/(3*PI)
      
     !Calculate interatomic spacing and electron relaxation time, R0/v_ave
     nion = nele/zbar
     R0 = (4*PI*nion/3)**(-1./3)
     relax_time = R0/sqrt(3*cond_boltz*xtemp/cond_mele)

     !Calculate Cowan melt temperature (and convert to K), LM84 Eq. 34
     !b = 0.6 * sim_zbarTarg**(1./9)        
     b = 0.6 * zbarTarg**(1./9)
     !xi = 9.0 * sim_zbarTarg**(0.3) * xden / sim_abarTarg        
     xi = 9.0 * zbarTarg**(0.3) * xden / abarTarg
     Tmelt = 0.32 * (xi/(1.+xi))**4 * xi**(2.*b-2./3) * 300 * 38.7

     if (xtemp > Tmelt) then
     tau = max(tau, relax_time)
     end if

     ! LM84, Eq. 23b
     cond = nele*cond_boltz*cond_boltz*xtemp*tau*A_beta/cond_mele

     vecLen = 1
     mode = MODE_DENS_TEMP
     eos_arr(EOS_TEMP) = xtemp
     eos_arr(EOS_DENS) = xden

     mask = .false.
     mask(EOS_CV)  = .true.
     mask(EOS_DET) = .true.

     call Eos(mode,vecLen,eos_arr,massfrac,mask)
     diff_coeff = cond/(xden*eos_arr(EOS_CV))
  
  else
     cond = 0.0
     diff_coeff = 0.0
  end if

end subroutine Conductivity
