!!****if* source/Simulation/SimulationMain/unitTest/Laser_quadraticTube/testI/sim_doAnalysis
!!
!!  NAME 
!!
!!   sim_doAnalysis
!!
!!  SYNOPSIS
!!
!!   sim_doAnalysis (logical (out) :: perfect)
!!
!!  DESCRIPTION
!!
!!   This routine does the analysis of the launched rays. Redirects according to the grid
!!   geometry present.
!!
!! ARGUMENTS
!!
!!  perfect : the success indicator for the simulation
!!
!!***

subroutine sim_doAnalysis (perfect)

  use Simulation_data,  ONLY : sim_geometry,      &
                               GRID_3DCARTESIAN,  &
                               GRID_2DCARTESIAN,  &
                               GRID_2DCYLINDRICAL

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

  logical, intent (out) :: perfect
!
!
!    ...Select the appropriate routine.
!
!
  select case (sim_geometry)

  case (GRID_2DCARTESIAN , GRID_2DCYLINDRICAL)

    call sim_doAnalysis2DRec (perfect)

  case (GRID_3DCARTESIAN)

    call sim_doAnalysis3DRec (perfect)

  case default

    call Driver_abortFlash ('[sim_doAnalysis] ERROR: Unknown geometry!')

  end select
!
!
!     ...Ready!
!
!
  return
end subroutine sim_doAnalysis
