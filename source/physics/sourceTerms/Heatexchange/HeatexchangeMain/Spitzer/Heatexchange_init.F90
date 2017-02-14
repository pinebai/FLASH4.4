!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Spitzer/Heatexchange_init
!!
!! NAME
!!
!!  Heatexchange_init
!!
!! SYNOPSIS
!!
!!  call Heatexchange_init(logical, intent(IN)  :: restart)
!!
!! DESCRIPTION
!!
!! Initialize Spitzer Heatexchange implementation. This just involves
!! reading a bunch of runtime parameters
!!
!! ARGUMENTS
!!
!!   restart : Indicates restart
!!
!!***

subroutine Heatexchange_init(restart)
  use Runtimeparameters_interface, ONLY: RuntimeParameters_get
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype
  use Heatexchange_data
  implicit none

#include "constants.h"
#include "Flash.h"

  
  logical, intent(IN) :: restart


  call Driver_getMype(MESH_COMM, hx_meshMe)
  hx_restart = restart

  call RuntimeParameters_get("useHeatexchange", hx_useHeatexchange)
  call RuntimeParameters_get("hx_logLevel"    , hx_logLevel)
  call RuntimeParameters_get("smalle", hx_smallE)

  call RuntimeParameters_get("hx_ieTimeCoef", hx_ieTimeCoef)
  call RuntimeParameters_get("hx_dtFactor", hx_dtFactor)
  call RuntimeParameters_get("hx_relTol", hx_relTol)
  if (hx_relTol .LT. 0.0) then
     call RuntimeParameters_get("eos_tolerance", hx_relTol)
  end if

  ! Set physical constants:
  call PhysicalConstants_get("electron mass",hx_mele)
  call PhysicalConstants_get("Boltzmann",hx_boltz)
  call PhysicalConstants_get("electron charge",hx_qele)
  call PhysicalConstants_get("Planck",hx_hbar)
  hx_hbar = hx_hbar/(2.0*PI)
  call PhysicalConstants_get("Avogadro", hx_navo)


end subroutine Heatexchange_init
