!!****if* source/Simulation/SimulationMain/unitTest/Multipole/sim_solveGridPoisson
!!
!!  NAME 
!!
!!   sim_solveGridPoisson
!!
!!  SYNOPSIS
!!
!!   sim_solveGridPoisson ()
!!
!!  DESCRIPTION
!!
!!   This routine computes the numerical gravitational potential for the Maclaurin
!!   problem. Supported boundary conditions are isolated only. The same boundary conditions
!!   are applied in all directions.
!!
!! ARGUMENTS
!!
!!***

subroutine sim_solveGridPoisson ()

  use Simulation_data,  ONLY : sim_pi,    &
                               sim_Newton

  use Timers_interface, ONLY : Timers_start,  &
                               Timers_stop

  use Grid_interface,   ONLY : GRID_PDE_BND_ISOLATED,  &
                               Grid_solvePoisson

  implicit none

#include "Flash.h"
#include "constants.h"
#include "Flash_mpi.h"

  integer :: ierr
  real    :: PoissonFactor

  integer :: bcTypes  (6)
  real    :: bcValues (2,6)
!
!
!   ...Preliminary chores.
!
!    
  call Timers_start ("MPI barrier in sim_solveGridPoisson")
  call MPI_Barrier  (MPI_COMM_WORLD, ierr)
  call Timers_stop  ("MPI barrier in sim_solveGridPoisson")

  bcTypes  = GRID_PDE_BND_ISOLATED       ! isolated in all 6 directions
  bcValues = 0.

  call Timers_start ("Solve grid Poisson Multipole")
!
!
!   ...Obtain the numerical gravitational solution.
!
!  
  PoissonFactor = 4. * sim_pi * sim_Newton

  call Grid_solvePoisson (GPOT_VAR,    &
                          DENS_VAR,    &
                          bcTypes,     &
                          bcValues,    &
                          PoissonFactor)

  call Timers_stop ("Solve grid Poisson Multipole")
!
!
!   ...Ready!
!
!    
  return
end subroutine sim_solveGridPoisson
