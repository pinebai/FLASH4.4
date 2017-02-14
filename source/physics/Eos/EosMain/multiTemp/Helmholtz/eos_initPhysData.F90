!!****if* source/physics/Eos/EosMain/multiTemp/Helmholtz/eos_initPhysData
!!
!! NAME
!!
!!  eos_initPhysData
!!
!! 
!! SYNOPSIS
!!
!!  call eos_initPhysData()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars used
!!  by the EOS implementation for physical constants,
!!  based on the runtime parameters and physical
!!  constants facilities. This is the version for a multiTemp
!!  Gamma law implementation with a single fluid.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!  
!!
!! NOTES
!!
!!
!!***
subroutine eos_initPhysData()

!!$  use Eos_data, ONLY: eos_gasConstant
!!$  use eos_helmConstData, ONLY: eos_avo, eos_avoinv, eos_kergavo, eos_ao3, &
!!$       eos_kerg, eos_sioncon, eos_h, eos_hbar, &
!!$       eos_c, eos_ssol, eos_asol, &
!!$       eos_amu
!!$!  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
!!$  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
!!$!  use Driver_interface, ONLY: Driver_abortFlash
!!$
 
  implicit none

#include "constants.h"
!#include "Flash.h"
!#include "Eos.h"  

!!$  real :: sb, c
!!$
!!$  !  call PhysicalConstants_get("ideal gas constant", eos_gasConstant)
!!$
!!$  !  call PhysicalConstants_get("electron mass",eos_eMass) !or value from eos_helmConstData?
!!$
!!$  !  call PhysicalConstants_get("electron mass",eos_eMassInUAmu,unitMass="amu")
!!$
!!$
!!$  !! DEV:  NOTE -- We should verify that unitsystem is "cgs" !!
!!$  call PhysicalConstants_get("Planck", eos_h)
!!$  eos_hbar    = 0.5 * eos_h / PI
!!$  call PhysicalConstants_get("Boltzmann", eos_kerg)
!!$
!!$  call PhysicalConstants_get("Avogadro", eos_avo)
!!$  eos_avoinv  = 1.0e0/eos_avo
!!$  eos_kergavo = eos_kerg * eos_avo
!!$  eos_amu = eos_avoinv          !!DEV: only valid in "cgs", not "MKS" units!
!!$
!!$  ! Compute radiation constant:
!!$  call PhysicalConstants_get("speed of light", c)
!!$  call PhysicalConstants_get("Stefan-Boltzmann", sb)
!!$  eos_ao3 = 4*sb/(3*c)
!!$
!!$  eos_sioncon = (2.0e0 * PI * eos_amu * eos_kerg)/(eos_h*eos_h)
!!$  return

end subroutine eos_initPhysData
