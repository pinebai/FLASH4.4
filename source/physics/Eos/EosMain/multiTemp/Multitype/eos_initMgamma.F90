!!****if* source/physics/Eos/EosMain/multiTemp/Multitype/eos_initMgamma
!!
!! NAME
!!
!!  eos_initMgamma
!!
!! SYNOPSIS
!!
!!  call eos_initMgamma()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars and arrays needed
!!  by the EOS unit from the runtime parameters and physical constants.
!!  This version is for use when multiple species
!!  are present. The gamma values for different species are obtained from
!!  the Multispecies unit, initialized in Simulation_initSpecies.F90 .
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos for multiple
!!   species with different abundances. 
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might over write these values with the 
!!   flash.par values for your specific run.  
!!
!!   eos_singleSpeciesA[Real]  -- Nucleon number for the gas, default 1.00
!!   eos_singleSpeciesZ[Real]  -- Proton number for the gas, default 1.00
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES and similar modes.
!!                                 If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!***

subroutine eos_initMgamma()

  use Eos_data, ONLY : eos_gasConstant, eos_smalle, eos_eintSwitch, eos_type, &
       eos_tol, eos_maxNewton, eos_largeT, eos_forceConstantInput, &
       eos_combinedTempRule, eos_entrEleScaleChoice
  use eos_mgammaData, ONLY : eos_gc, eos_gammam1j, eos_ggprodj, eos_ggprodinvj, eos_gam1invj, &
       eos_eMass, eos_eMassInUAmu, &
       eos_gammaEle, eos_gammam1Ele, eos_gammam1Rad, &
       eos_maxFactorUp, eos_maxFactorDown
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Multispecies_interface, ONLY : Multispecies_getProperty
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use eos_mtInterface, ONLY: eos_initPhysData

 
  implicit none

#include "Flash.h"
#include "Eos.h"
#include "Multispecies.h"


  integer :: ispecies

  call RuntimeParameters_get("gammaEle", eos_gammaEle)
  call RuntimeParameters_get('eos_largeT', eos_largeT)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get('eos_maxFactorUp', eos_maxFactorUp)
  call RuntimeParameters_get('eos_maxFactorDown', eos_maxFactorDown)

#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abortFlash("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif

#ifdef OLD_EINTNSWITCH
  call RuntimeParameters_get("eint1Switch",eos_eint1Switch)  
  call RuntimeParameters_get("eint2Switch",eos_eint2Switch)  
  call RuntimeParameters_get("eint3Switch",eos_eint3Switch)  
  if (eos_eint1Switch==-1.0) eos_eint1Switch = eos_eintSwitch
  if (eos_eint2Switch==-1.0) eos_eint2Switch = eos_eintSwitch
  if (eos_eint3Switch==-1.0) eos_eint3Switch = eos_eintSwitch
#endif
  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)
  call RuntimeParameters_get("eos_combinedTempRule",eos_combinedTempRule)
  call RuntimeParameters_get("eos_entrEleScaleChoice", eos_entrEleScaleChoice)
  call PhysicalConstants_get("electron mass",eos_eMass) !or value from eos_helmConstData?

  call PhysicalConstants_get("electron mass",eos_eMassInUAmu,unitMass="amu")

  ! Compute radiation constant, maybe other physical constants:
  call eos_initPhysData()

  eos_gammam1Ele = 1.0/(eos_gammaEle-1.0)
  eos_gammam1Rad = 3.0
  

  do ispecies = 1, NSPECIES
     call Multispecies_getProperty(SPECIES_BEGIN + ispecies - 1, GAMMA,  eos_gc(ispecies))
  end do

  ! Note that these are all ARRAYS of size NSPECIES
  eos_gammam1j(:)   = 1. / (eos_gc(:) - 1.)
  eos_ggprodj(:)    = eos_gammam1j(:) * eos_gasConstant
  eos_ggprodinvj(:) = 1. / eos_ggprodj(:)
  eos_gam1invj(:)   = 1. / eos_gammam1j(:)

  eos_type = EOS_MTYPE

  return
end subroutine eos_initMgamma
