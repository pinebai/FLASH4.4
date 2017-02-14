!!****if* source/Simulation/SimulationMain/Plasma/Simulation_init
!!
!! NAME
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!  Initializes all the parameters needed for simulating an ion beam
!!  into plasma. The simulation uses a Hyrbid-pic method for the computation
!!  
!! ARGUMENTS
!!  
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call RuntimeParameters_get('pt_picPmass_1', sim_pmass_1)
  call RuntimeParameters_get('pt_picPvelx_1', sim_pvelx_1)
  call RuntimeParameters_get('pt_picPvely_1', sim_pvely_1)
  call RuntimeParameters_get('pt_picPvelz_1', sim_pvelz_1)

  call RuntimeParameters_get('pt_picPmass_2', sim_pmass_2)

  call RuntimeParameters_get('sim_bx', sim_bx)
  call RuntimeParameters_get('sim_by', sim_by)
  call RuntimeParameters_get('sim_bz', sim_bz)

  sim_usw = (/ sim_pvelx_1, sim_pvely_1, sim_pvelz_1 /)
  sim_bsw = (/ sim_bx, sim_by, sim_bz /)
  ! SW electric field from E = -uxB in SW
  sim_esw(1) = sim_usw(2)*sim_bsw(3)-sim_usw(3)*sim_bsw(2)
  sim_esw(2) = sim_usw(3)*sim_bsw(1)-sim_usw(1)*sim_bsw(3)
  sim_esw(3) = sim_usw(1)*sim_bsw(2)-sim_usw(2)*sim_bsw(1)
  sim_esw = -sim_esw

end subroutine Simulation_init
