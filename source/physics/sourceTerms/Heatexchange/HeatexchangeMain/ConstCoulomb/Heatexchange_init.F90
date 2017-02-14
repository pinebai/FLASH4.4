!!****if* source/physics/sourceTerms/Heatexchange/HeatexchangeMain/ConstCoulomb/Heatexchange_init
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
!!
!! ARGUMENTS
!!
!!   restart : indicates restart
!!
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

  call RuntimeParameters_get("hx_coulombLog", hx_coulombLog)
  call RuntimeParameters_get("hx_couplingConst13", hx_c13)
  call RuntimeParameters_get("hx_couplingConst23", hx_c23)

  
  call RuntimeParameters_get("hx_dtFactor", hx_dtFactor)
  call RuntimeParameters_get("hx_relTol", hx_relTol)
  if (hx_relTol .LT. 0.0) then
     call RuntimeParameters_get("eos_tolerance", hx_relTol)
  end if

  call RuntimeParameters_get("eos_singleSpeciesA", hx_singleSpeciesA)
  call RuntimeParameters_get("eos_singleSpeciesZ", hx_singleSpeciesZ)

  call PhysicalConstants_get("Avogadro", hx_Avogadro)
  call PhysicalConstants_get("Boltzmann", hx_kBoltzmann)
  call PhysicalConstants_get("electron charge",hx_eleCharge)
  call PhysicalConstants_get("electron mass",hx_eMassInUAmu,unitMass="amu")


#if(0)
  Lamx = 1.5 * ne * kBoltzmann * diffT * nuj
#endif

end subroutine Heatexchange_init
