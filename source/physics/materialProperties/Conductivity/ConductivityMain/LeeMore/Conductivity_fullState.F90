!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/LeeMore/Conductivity_fullState
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
!! Computes the Lee & More (Phys. Fluids 27, 1273 1984) electron thermal conductivity 
!! for all materials. Hereafter this paper is referred to as LM84.
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
#include "Multispecies.h"

subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo, cond_hbar
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar, Eos_getTempData
  use Multispecies_interface, ONLY: Multispecies_getSumFrac, Multispecies_getSumInv

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

  real :: xtemp, xden, tion
  real :: nele, ll
  real :: abar, zbar, znucbar
  real :: zbarFrac

  real :: tau, A_beta, nion, relax_time, R0
  real :: mu_div_kT, Fermi_onehalf, n3D, y
  real :: b, xi, Tmelt

  real, parameter :: a1 = 13.5
  real, parameter :: a2 = 0.976
  real, parameter :: a3 = 0.437
  real, parameter :: b2 = 0.51
  real, parameter :: b3 = 0.126

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0
  
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
     call Driver_abortFlash("[Conductivity_fullState] LeeMore conductivity computes only electron conduction, not ion conduction")
#endif
  case(2)
#ifdef TELE_VAR
     tempToUse = TELE_VAR
#endif
  case(3)
#ifdef TRAD_VAR
     call Driver_abortFlash("[Conductivity_fullState] LeeMore conductivity does not work with radiation")
#endif
  end select

    if (cond_useConductivity) then
       call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)

       massfrac => solnVec(SPECIES_BEGIN:SPECIES_END)
       call Multispecies_getSumFrac(Z,zbarFrac,massfrac)
       znucbar = abar*zbarFrac   ! mean atomic number in cell

       xden = solnVec(DENS_VAR)
       xtemp = solnVec(tempToUse)
       tion = solnVec(TION_VAR)
       nele = zbar * xden * cond_navo / abar
       call cond_logLambda(xtemp, nele, tion, zbar, ll)

       call fermi_chem_potential(nele, xtemp, mu_div_kT)

       n3D = 2 * (cond_mele * cond_boltz * xtemp / (2 * PI * cond_hbar**2 ))**(3./2)
       Fermi_onehalf = nele/n3D

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
       b = 0.6 * znucbar**(1./9)       
       xi = 9.0 * znucbar**(0.3) * xden / abar
       Tmelt = 0.32 * (xi/(1.+xi))**4 * xi**(2*b-2./3) * 300 * 38.7

       if (xtemp > Tmelt) then
       tau = max(tau, relax_time)
       end if

       ! LM84, Eq. 23b
       isochoricCondLoc = nele*cond_boltz*cond_boltz*xtemp*tau*A_beta/cond_mele

       if (present(diffCoeff)) then

          eos_arr(EOS_DENS) = xden
          eos_arr(EOS_TEMPELE) = solnVec(TELE_VAR)
          eos_arr(EOS_TEMPION) = solnVec(TION_VAR)
          eos_arr(EOS_TEMPRAD) = solnVec(TRAD_VAR)

!          call Eos_getTempData(solnVec,eos_arr,MODE_DENS_TEMP_GATHER)

          mask = .false.
          mask(EOS_CVELE)  = .true.
          mask(EOS_DET) = .true.

          if (NSPECIES > 0) then
             !No need to calculate massfrac again, this was called earlier to get znucbar
             !massfrac => solnVec(SPECIES_BEGIN:SPECIES_END) 
             call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,massfrac,mask)
          else
             call Eos(MODE_DENS_TEMP_GATHER,1,eos_arr,mask=mask)
          end if
          diffCoeffLoc = isochoricCondLoc/(xden*eos_arr(EOS_CVELE))
       end if

    end if

#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc
  
end subroutine Conductivity_fullState
