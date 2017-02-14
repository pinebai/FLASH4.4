!!****if* source/Simulation/SimulationMain/StirFromFile/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!! DESCRIPTION
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!
!!***

subroutine Simulation_init()

  use Driver_interface, ONLY: Driver_getMype
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "Flash.h"
#include "constants.h"

  integer :: dr_myPE

  call RuntimeParameters_get('rho_ambient', sim_rhoAmbient)
  call RuntimeParameters_get('c_ambient', sim_cAmbient)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('magnetic', sim_magnetic)
  if (sim_magnetic) then
    call RuntimeParameters_get('MagField_z', sim_MagField_z)
  endif

  call Driver_getMype(GLOBAL_COMM, dr_myPE)

  if (dr_myPE .eq. MASTER_PE) then
    write(*,'(A)') 'Initializing the StirringFromFile turbulence setup...'
    if (sim_magnetic) then
#ifdef MAGZ_VAR
      print *, 'Constant magnetic field was set to ', sim_MagField_z, ' in z-direction.'
#endif
    endif
  endif

end subroutine Simulation_init
