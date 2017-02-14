!!****if* source/Simulation/SimulationMain/unitTest/PFFT_BlktriFD/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! ARGUMENTS
!!
!!   
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for INS-isotropic turbulence problem.
!!
!!***

subroutine Simulation_init()

  use Simulation_data

  use Driver_interface, ONLY : Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  implicit none

#include "constants.h"
#include "Flash.h"

  

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get('xmin',    sim_xMin)
  call RuntimeParameters_get('ymin',    sim_yMin)
  call RuntimeParameters_get('zmin',    sim_zMin)
  call RuntimeParameters_get('xmax',    sim_xMax)
  call RuntimeParameters_get('ymax',    sim_yMax)
  call RuntimeParameters_get('zmax',    sim_zMax)

  call RuntimeParameters_get("waven_x",  pfb_waven_x)
  call RuntimeParameters_get("waven_y",  pfb_waven_y)
  call RuntimeParameters_get("waven_z",  pfb_waven_z)
  call RuntimeParameters_get("alpha_x",  pfb_alpha_x)

  if (sim_meshMe .eq. 0 ) then
     write(*,*) 'waven_x =',pfb_waven_x
     write(*,*) 'waven_y =',pfb_waven_y
     write(*,*) 'waven_z =',pfb_waven_z
     write(*,*) 'alpha_x =',pfb_alpha_x
  endif

  sim_gCell = .true.

end subroutine Simulation_init
