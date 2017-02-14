!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/LeeMore/fermi_chem_potential
!!
!! NAME
!!
!!  fermi_chem_potential
!!
!! SYNOPSIS
!!
!!  call fermi_chem_potential(     real(IN)  :: nele,
!!                                 real(IN)  :: tele,
!!                                 real(OUT) :: mu_div_kT)
!!
!! DESCRIPTION
!!
!! Computes the chemical potential divided by the thermal energy for a fermi gas
!! by a handy fitting function to the inverse of the Fermi(j = 1/2) integral 
!! which appears in Zimmermann, R., Many-Particle Theory of Highly Excited 
!! Semiconductors, Teubner Verlagsgesellschaft, Leipzig (1988) p. 150. 
!! 
!!
!! ARGUMENTS
!!
!!   nele      :   electron number density 
!!   tele      :   electron temperature 
!!   mu_div_kT :   chemical potential divided by k T
!!
!! SEE ALSO
!!
!!  Conductivity 
!!***

#include "constants.h"

  subroutine fermi_chem_potential(nele, tele, mu_div_kT)
    use Conductivity_data, ONLY: cond_mele, cond_boltz, cond_hbar, cond_qele
    implicit none

    real, intent(in)  :: nele ! electron number density [cm^-3]
    real, intent(in)  :: tele ! electron temperature [K]
    real, intent(out) :: mu_div_kT   ! chemical potential divided by k T [unitless]

    real :: y, n3D

    n3D = 2 * (cond_mele * cond_boltz * tele / (2*PI*cond_hbar**2) )**(3./2)

    y = nele / n3D

    if (y < 5.5) then
    mu_div_kT = log(y) + 0.3536*y - 0.00495*y**2 + 0.000125*y**3
    else if (y > 5.5) then
    mu_div_kT = 1.209*y**(2./3) - 0.6803*y**(-2./3) - 0.85/(y*y)
    end if

  end subroutine fermi_chem_potential

