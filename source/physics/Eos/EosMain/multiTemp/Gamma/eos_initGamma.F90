!!****if* source/physics/Eos/EosMain/multiTemp/Gamma/eos_initGamma
!!
!! NAME
!!
!!  eos_initGamma
!!
!! SYNOPSIS
!!  
!!  call  eos_initGamma()
!!                 
!!
!! DESCRIPTION
!!
!!  This routine does ideal gamma law specific initialization
!!
!!  ARGUMENTS
!!
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos. 
!!   Particular implementations (Gamma,Helmholz,etc) of the unit
!!   define their own runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might over write these values with the 
!!   flash.par values for your specific run.  
!!
!!   gamma[Real]   -- Ratio of specific heats, default 1.6667
!!   eos_singleSpeciesA[Real]  -- Nucleon number for the gas, default 1.00
!!   eos_singleSpeciesZ[Real]  -- Proton number for the gas, default 1.00
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES. If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!
!!  NOTES
!!
!!  Gamma law Eos defines two mesh-based parameters GAMC_VAR and GAME_VAR in Flash.h
!!
!!***

#include "Eos.h"
subroutine eos_initGamma()

  use Eos_data, ONLY : eos_gasConstant, eos_eintSwitch, eos_type, &
       eos_tol, eos_maxNewton, eos_forceConstantInput, &
       eos_combinedTempRule, eos_entrEleScaleChoice
  use eos_idealGammaData, ONLY : eos_gamma, eos_gammam1, &
       eos_singleSpeciesA, eos_singleSpeciesZ, &
       eos_eMass, eos_eMassInUAmu, &
       eos_gammaIon, eos_gammaEle, eos_gammaRad, eos_gammam1Ion, eos_gammam1Ele, eos_gammam1Rad

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use eos_mtInterface, ONLY: eos_initPhysData

  implicit none

#include "Flash.h"
#include "Eos.h"  

  eos_type=EOS_GAM

  call PhysicalConstants_get("ideal gas constant", eos_gasConstant)

  call RuntimeParameters_get("gamma", eos_gamma)
  call RuntimeParameters_get("gammaIon", eos_gammaIon)
  call RuntimeParameters_get("gammaEle", eos_gammaEle)
  call RuntimeParameters_get("gammaRad", eos_gammaRad)
  call RuntimeParameters_get("eos_singleSpeciesA", eos_singleSpeciesA)
  call RuntimeParameters_get("eos_singleSpeciesZ", eos_singleSpeciesZ)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
#ifndef EINT_VAR
  if (eos_eintSwitch > 0.0) then
     call Driver_abortFlash("[Eos_init] eintSwitch is nonzero, but EINT_VAR not defined!")
  end if
#endif
  call RuntimeParameters_get("eos_forceConstantInput",eos_forceConstantInput)
  call RuntimeParameters_get("eos_combinedTempRule",eos_combinedTempRule)
  call RuntimeParameters_get("eos_entrEleScaleChoice", eos_entrEleScaleChoice)
  call PhysicalConstants_get("electron mass",eos_eMass) !or value from eos_helmConstData?

  call PhysicalConstants_get("electron mass",eos_eMassInUAmu,unitMass="amu")

  ! Compute radiation constant, maybe other physical constants:
  call eos_initPhysData()

  eos_gammam1 = 1.0/(eos_gamma-1.0)
  eos_gammam1Ion = 1.0/(eos_gammaIon-1.0)
  eos_gammam1Ele = 1.0/(eos_gammaEle-1.0)
  eos_gammam1Rad = 1.0/(eos_gammaRad-1.0)
  

  return
end subroutine eos_initGamma
