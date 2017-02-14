!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_launchRays
!!
!! NAME 
!!
!!  sim_launchRays
!!
!! SYNOPSIS
!!
!!  sim_launchRays ()
!!
!! DESCRIPTION
!!
!!  This routine launches the rays. Redirects according to the grid geometry present.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_launchRays ()

  use Simulation_data,  ONLY : sim_geometry,      &
                               GRID_3DCARTESIAN,  &
                               GRID_2DCARTESIAN,  &
                               GRID_2DCYLINDRICAL

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
!
!
!    ...Select the appropriate routine.
!
!
  select case (sim_geometry)

  case (GRID_2DCARTESIAN , GRID_2DCYLINDRICAL)

    call sim_launchRays2DRec ()

  case (GRID_3DCARTESIAN)

    call sim_launchRays3DRec ()

  case default

    call Driver_abortFlash ('[sim_launchRays] ERROR: Unknown geometry!')

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine sim_launchRays
