!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock (integer (in) :: blockID)
!!
!! DESCRIPTION
!!
!!  Initializes the data for a specified block needed to do the laser quadratic tube unit test.
!!  Redirects according to the grid geometry present.
!!
!! ARGUMENTS
!!
!!  blockID : The block ID number to be initialized
!!
!!***

subroutine Simulation_initBlock (blockID)

  use Simulation_data,   ONLY : sim_geometry,      &
                                GRID_3DCARTESIAN,  &
                                GRID_2DCARTESIAN,  &
                                GRID_2DCYLINDRICAL

  use Driver_interface,  ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent (in) :: blockID
!
!
!    ...Select the appropriate routine.
!
!
  select case (sim_geometry)

  case (GRID_2DCARTESIAN , GRID_2DCYLINDRICAL)

    call sim_initBlock2DRec (blockID)

  case (GRID_3DCARTESIAN)

    call sim_initBlock3DRec (blockID)

  case default

    call Driver_abortFlash ('Laser quadratic tube: unsupported geometry (Simulation_initBlock)!')

  end select
!
!
!    ...Ready!
!
!
  return
end subroutine Simulation_initBlock
