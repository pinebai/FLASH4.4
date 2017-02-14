!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Constant/Heatexchange_init
!!
!! NAME
!!
!!  Heatexchange_init
!!
!! SYNOPSIS
!!
!!  call Heatexchange_init(logical(IN)  :: restart)
!!
!! DESCRIPTION
!!
!!   Initializes the Heatexchange unit.
!!
!! ARGUMENTS
!!
!!   restart : indicates restart
!!
!!***

subroutine Heatexchange_init(restart)
  use Runtimeparameters_interface, ONLY: RuntimeParameters_get
  use PhysicalConstants_interface, ONLY: PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype
  use Heatexchange_data
  implicit none
  
#include "constants.h"

  logical, intent(IN) :: restart

  call Driver_getMype(MESH_COMM, hx_meshMe)
  hx_restart = restart

  call RuntimeParameters_get("useHeatexchange", hx_useHeatexchange)
  call RuntimeParameters_get("smalle", hx_smallE)

  call RuntimeParameters_get("hx_couplingConst12", hx_c12)
  call RuntimeParameters_get("hx_couplingConst13", hx_c13)
  if (hx_c13 == -1) then
     hx_c13 = hx_c12
  end if
  call RuntimeParameters_get("hx_couplingConst23", hx_c23)
  if (hx_c23 == -1) then
     hx_c23 = hx_c12
  end if

  
  call RuntimeParameters_get("hx_dtFactor", hx_dtFactor)
  call RuntimeParameters_get("hx_relTol", hx_relTol)
  if (hx_relTol .LT. 0.0) then
     call RuntimeParameters_get("eos_tolerance", hx_relTol)
  end if

  ! Set some physical constants:
  call PhysicalConstants_get("Boltzmann", hx_kBoltzmann)

  call RuntimeParameters_get("useRadTrans", hx_useRadTrans)

end subroutine Heatexchange_init
