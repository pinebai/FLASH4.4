!!****if* source/physics/Eos/EosMain/multiTemp/Eos_init
!!
!! NAME
!!
!!  Eos_init
!!
!! 
!! SYNOPSIS
!!
!!  call Eos_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes various scalars needed
!!  by the EOS unit from the runtime parameters and physical
!!  constants facilities. This is the version for simple Gamma law
!!  implementation with a single fluid.
!!
!! ARGUMENTS
!!
!!  none
!!
!! PARAMETERS
!!  
!!   These are the runtime parameters used in Gamma law Eos. 
!!   Particular implementations (Gamma,Helmholz,etc) of the unit
!!   define their own runtime parameters.
!!
!!   To see the default parameter values and all the runtime parameters
!!   specific to your simulation check the "setup_params" file in your
!!   object directory. You might overwrite these values with the 
!!   flash.par values for your specific run.  
!!
!!   gamma[Real]   -- Ratio of specific heats, default 1.6667
!!   eos_singleSpeciesA[Real]  -- Nucleon number for the gas, default 1.00
!!   eos_singleSpeciesZ[Real]  -- Proton number for the gas, default 1.00
!!   eos_logLevel[Integer]     -- Verbosity of some warnings
!!   eos_tolerance[Real, 1.0e-8] -- Convergence tolerance for the Newton-Rhapson
!!               iterations
!!   eos_maxNewton[Integer, 50] -- Maximum number of Newton-Raphson iterations
!!   eos_forceConstantInput     -- This switch forces the Eos implementation
!!                                 to never modify EINT or PRES in MODE_DENS_EI
!!                                 and MODE_DENS_PRES. If this is .false. (the
!!                                 default), calls to Eos may slightly modify
!!                                 these input variables in order to preserve
!!                                 thermodynamic equilibrium.
!!   eos_smallEion
!!   eos_smallEele
!!   eos_smallErad
!!
!!  NOTES
!!
!!   Eos defines at least two mesh-based variables GAMC_VAR and GAME_VAR in Flash.h
!!
!!***
subroutine Eos_init()

  use Eos_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Driver_interface, ONLY: Driver_abortFlash
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs
  use eos_mtInterface, ONLY: eos_initPhysData
  use eos_localInterface, ONLY : eos_initMgamma, eos_initHelmholtz,&
       eos_initMtemp,eos_initTabulated, eos_initGamma
  use Timers_interface, ONLY : Timers_start, Timers_stop  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"  

  call Timers_start("Eos_init")

  ! Everybody should know this
  call Driver_getMype(MESH_COMM,eos_meshMe)
  call Driver_getMype(GLOBAL_COMM, eos_globalMe)
  call Driver_getNumProcs(MESH_COMM,eos_meshNumProcs)

  call PhysicalConstants_get("ideal gas constant", eos_gasConstant)

  call RuntimeParameters_get("eos_singleSpeciesA", eos_singleSpeciesA)
  call RuntimeParameters_get("eos_singleSpeciesZ", eos_singleSpeciesZ)
  call RuntimeParameters_get("eos_logLevel", eos_logLevel)
  call RuntimeParameters_get('smallt', eos_smallt)
  call RuntimeParameters_get('smlrho', eos_smallRho)
  call RuntimeParameters_get('eos_tolerance', eos_tol)
  call RuntimeParameters_get('eos_maxNewton', eos_maxNewton)
  call RuntimeParameters_get("smalle",eos_smalle)
  call RuntimeParameters_get("eintSwitch",eos_eintSwitch)
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
  call RuntimeParameters_get("eos_smallEion", eos_smallEion)
  call RuntimeParameters_get("eos_smallEele", eos_smallEele)
  call RuntimeParameters_get("eos_smallErad", eos_smallErad)

  ! Compute radiation constant, maybe other physical constants:
  call eos_initPhysData()


  call eos_fillMapLookup()

  call eos_initGamma()
  call eos_initMgamma()
  call eos_initHelmholtz()
  call eos_initMtemp()
  call eos_initTabulated()

  call Timers_stop("Eos_init")

  return

end subroutine Eos_init
