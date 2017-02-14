!!***if* source/physics/Cosmology/CosmologyMain/Cosmology_init
!!
!! NAME
!!
!!  Cosmology_init
!!
!!
!! SYNOPSIS
!!
!!  Cosmology_init(logical(IN) :: restart)
!!  
!!
!! DESCRIPTION
!! 
!!  Initialize unit scope variables which are typically the runtime parameters.
!!  This must be called once by Driver_initFlash.F90 first. Calling multiple
!!  times will not cause any harm but is unnecessary.
!!
!! ARGUMENTS
!!
!!  
!!  restart -- note if we are restarting. 
!!
!! NOTE
!!
!!  From a "cold start" this unit has to determine the cosmological scaling 
!!  and will solve the Friedmann equations to do so.  On a restart, the
!!  correct values will be read from the checkpoint file
!!
!!***
Subroutine Cosmology_init( restart)

  use Cosmology_data

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use IO_interface, ONLY : IO_getScalar
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs  
  implicit none

#include "Flash.h"
#include "constants.h"

  real :: dt_init
  real :: redshift, redshiftOld
  logical, intent(IN) :: restart
 
  real :: initRedshift
  real :: simTime

  call Driver_getMype(MESH_COMM,csm_meshMe)
  call Driver_getNumProcs(MESH_COMM,csm_meshNumProcs)

  call RuntimeParameters_get("smlrho", csm_smlrho)
  call RuntimeParameters_get("smalle", csm_smalle)
  call RuntimeParameters_get("smallp", csm_smallp)
#ifdef EINT_VAR
  call RuntimeParameters_get("eintSwitch",csm_eintSwitch)
#endif
  call RuntimeParameters_get("HubbleConstant", csm_hubble)
  call RuntimeParameters_get("OmegaMatter", csm_omega)
  call RuntimeParameters_get("OmegaBaryon", csm_baryon)
  call RuntimeParameters_get("CosmologicalConstant", csm_lambda)
  call RuntimeParameters_get("MaxScaleChange", csm_maxScaleChange)
  call RuntimeParameters_get("zFinal", csm_redshiftFinal)
  
  call PhysicalConstants_get("Newton",csm_newton)
  call PhysicalConstants_get("speed of light",csm_c)
  
  call RuntimeParameters_get("dtinit", dt_init)

  call RuntimeParameters_get("computeRedshiftOnly", csm_computeRedshiftOnly)

  csm_curv = csm_omega + csm_lambda - 1.

  if(.NOT. restart) then
     
     call RuntimeParameters_get("zInitial", initRedshift)
     csm_scaleFactor = 1./(1.+ initRedshift)
     csm_oldScaleFactor = csm_scaleFactor * (1.-csm_maxScaleChange)
     
  else 
     
     ! If we are restarting, solving the Friedmann equation is redunant

     call IO_getScalar('scaleFactor', csm_scaleFactor)
     call IO_getScalar('oldScaleFactor',csm_oldScaleFactor)
     
     !Time is read in by Driver on restart
     !Friedmann equation already "solved"

  end if

end Subroutine Cosmology_init
