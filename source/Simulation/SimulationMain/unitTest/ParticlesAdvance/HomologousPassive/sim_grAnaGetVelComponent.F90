!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/sim_grAnaGetVelComponent
!!
!!  SYNOPSIS
!!
!!  sim_grAnaGetVelComponent( real(out)   :: result,
!!                             integer(in) :: varGrid,
!!                             integer(in) :: xCoord,
!!                             integer(in) :: yCoord,
!!                             integer(in) :: zCoord,
!!                             real(in)    :: t)
!!
!!
!! DESCRIPTION
!!  function to compute the value of velocity component analytically from 
!!  the grid
!!
!! ARGUMENTS
!!
!! result   : argument that returns the computed result
!! varGrid  : the index into the grid data structure
!! xCoord   : particle property with x component of velocity
!! yCoord   : particle property with y component of velocity
!! zCoord   : particle property with z component of velocity
!! t        : time
!!
!!***
subroutine sim_grAnaGetVelComponent(result,&
     varGrid,xCoord,yCoord,zCoord,t)

  use Driver_interface, ONLY : Driver_abortFlash
  use Simulation_data, ONLY : sim_a0, sim_a1
  implicit none

#include "Flash.h"

  real, INTENT(out)           :: result
  integer, INTENT(in)           :: varGrid
  real, INTENT(in)           :: xCoord, yCoord, zCoord, t

  real :: expansionFactor, relevantCoord

  if (varGrid == VELX_VAR) then
     relevantCoord = xCoord
  else if (varGrid == VELY_VAR) then
     relevantCoord = yCoord
  else if (varGrid == VELZ_VAR) then
     relevantCoord = zCoord
  else
     print*, "sim_grAnaGetVelComponent should not have been called for this varGrid:", varGrid
     call Driver_abortFlash("sim_grAnaGetVelComponent should not have been called for this varGrid!")
  end if

  expansionFactor = sim_a0 + sim_a1 * t
  result = relevantCoord * expansionFactor

end subroutine sim_grAnaGetVelComponent
