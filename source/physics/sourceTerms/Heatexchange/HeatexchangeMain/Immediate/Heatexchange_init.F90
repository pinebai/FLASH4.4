!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/Immediate/Heatexchange_init
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
!!   restart : indicates whether the run is restarting from a checkpoint
!!             (ignored in this implementation)
!!
!! PARAMETERS
!!
!!  useHeatexchange -- Boolean, True.  Turns on Heatexchange unit
!!  hx_applyToRadiation -- Boolean, False.  Does the Heatexchange unit handle radiation?
!!
!!***

subroutine Heatexchange_init(restart)
  use Runtimeparameters_interface, ONLY: RuntimeParameters_get
  use Driver_interface, ONLY : Driver_getMype
  use Heatexchange_data, ONLY : hx_useHeatexchange, hx_applyToRadiation, &
                                hx_globalMe, hx_restart
  implicit none
  
#include "constants.h"

  logical, intent(IN) :: restart

  call Driver_getMype(GLOBAL_COMM, hx_globalMe)
  hx_restart = restart

  call RuntimeParameters_get("useHeatexchange", hx_useHeatexchange)
  call RuntimeParameters_get("hx_applyToRadiation", hx_applyToRadiation)

  if (hx_useHeatexchange) then
     if (hx_globalMe == MASTER_PE) then
        print *,"NOTE: Using Immediate Heatexchange implementation."
     end if
  end if

end subroutine Heatexchange_init
