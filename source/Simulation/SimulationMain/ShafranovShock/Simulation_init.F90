!!****if* source/Simulation/SimulationMain/ShafranovShock/Simulation_init
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
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the Sod shock tube
!!  problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_rhoLeft    Density in the left part of the grid
!!  sim_rhoRight   Density in the right part of the grid
!!  sim_pLeft      Pressure  in the left part of the grid
!!  sim_pRight     Pressure  in the righ part of the grid
!!  sim_uLeft      fluid velocity in the left part of the grid
!!  sim_uRight     fluid velocity in the right part of the grid
!!  sim_xangle     Angle made by diaphragm normal w/x-axis (deg)
!!  sim_ yangle    Angle made by diaphragm normal w/y-axis (deg)
!!  sim_posnR      Point of intersection between the shock plane and the x-axis
!!
!!***

subroutine Simulation_init()
  
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Driver_interface,  ONLY : Driver_abortFlash, Driver_getMype
  implicit none
#include "Flash.h"
#include "constants.h"
  

  integer :: i

  if (NDIM /= 1) then
     call Driver_abortFlash('[Simulation_init] ERROR This problem is defined in 1D only!')
  end if

  call RuntimeParameters_get('smallx', sim_smallX)

  call RuntimeParameters_get('sim_DataPoints', sim_DataPoints)
  call RuntimeParameters_get('sim_InitData'  , sim_InitData)  
  call RuntimeParameters_get('sim_ShockSpeed', sim_ShockSpeed)

  call Logfile_stamp( "initializing Shafranov problem",  &
       "[Simulation_init]")

  allocate (sim_x(sim_DataPoints))
  allocate (sim_velx(sim_DataPoints))
  allocate (sim_tele(sim_DataPoints))
  allocate (sim_rho(sim_DataPoints))
  allocate (sim_tion(sim_DataPoints))
  allocate (sim_pele(sim_DataPoints))
  allocate (sim_pion(sim_DataPoints))
  allocate (sim_sele(sim_DataPoints))

  open (unit = 4, file=sim_InitData)
  do i=1, sim_DataPoints
      read (4,*) sim_x(i), sim_velx(i), sim_tele(i), sim_rho(i), &
                 sim_tion(i), sim_pele(i), sim_pion(i), sim_sele(i)
  enddo 
  close(4)

  call RuntimeParameters_get('sim_maxTol', sim_maxTol)
  call Driver_getMype(MESH_COMM, sim_meshMe)
     
end subroutine Simulation_init







