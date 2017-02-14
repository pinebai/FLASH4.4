!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/sim_ptAnaGetVelComponent
!!
!! NAME
!!
!!  sim_ptAnaGetVelComponent
!!
!! SYNOPSIS
!!
!!  call sim_ptAnaGetVelComponent(real, INTENT(inout),dimension(NPART_PROPS,maxParticlesPerProc)  :: particles,
!!                                integer, INTENT(in)  :: maxparticlesperproc,
!!                                integer, INTENT(in)  :: vargrid,
!!                                integer, INTENT(in)  :: xpart,
!!                                integer, INTENT(in)  :: ypart,
!!                                integer, INTENT(in)  :: zpart,
!!                                real, INTENT(in)  :: t,
!!                                integer, INTENT(in)  :: proppart,
!!                                integer, INTENT(in)  :: k)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   particles : 
!!
!!   maxparticlesperproc : 
!!
!!   vargrid : 
!!
!!   xpart : 
!!
!!   ypart : 
!!
!!   zpart : 
!!
!!   t : 
!!
!!   proppart : 
!!
!!   k : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine sim_ptAnaGetVelComponent(particles,maxParticlesPerProc,&
     varGrid,xPart,yPart,zPart,t,propPart,k)

  use Driver_interface, ONLY : Driver_abortFlash
  use Simulation_data, ONLY : sim_a0, sim_a1
  implicit none

#include "Flash.h"

  integer, INTENT(in) :: maxParticlesPerProc
  real, INTENT(inout),dimension(NPART_PROPS,maxParticlesPerProc) :: particles
  integer, INTENT(in)           :: varGrid, propPart
  integer, INTENT(in)           :: xPart, yPart, zPart, k
  real, INTENT(in)           :: t

  real :: expansionFactor, relevantCoord

  if (varGrid == VELX_VAR) then
     relevantCoord = particles(xPart,k)
  else if (varGrid == VELY_VAR) then
     relevantCoord = particles(yPart,k)
  else if (varGrid == VELZ_VAR) then
     relevantCoord = particles(zPart,k)
  else
     print*, "sim_ptAnaGetVelComponent should not have been called for this varGrid:", varGrid
     call Driver_abortFlash("sim_ptAnaGetVelComponent should not have been called for this varGrid!")
  end if

  expansionFactor = sim_a0 + sim_a1 * t
  particles(propPart,k) = relevantCoord * expansionFactor

end subroutine sim_ptAnaGetVelComponent
