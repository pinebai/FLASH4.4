!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/LeeMore/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  Conductivity_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***
subroutine Conductivity_init
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_hbar, cond_qele, cond_navo
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "Flash.h"
#include "constants.h"

  

  call RuntimeParameters_get("useConductivity", cond_useConductivity)

  ! Set physical constants:
  call PhysicalConstants_get("electron mass",cond_mele)
  call PhysicalConstants_get("Boltzmann",cond_boltz)
  call PhysicalConstants_get("electron charge",cond_qele)
  call PhysicalConstants_get("Planck",cond_hbar)
  cond_hbar = cond_hbar/(2.0*PI)
  call PhysicalConstants_get("Avogadro", cond_navo)

end subroutine Conductivity_init

