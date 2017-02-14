!!****ih* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/sim_interface
!!
!! NAME
!!
!!  sim_interface
!!
!! SYNOPSIS
!!
!!  use sim_interface
!!
!! DESCRIPTION
!!
!!  Interfaces for some subprograms private to the unit test
!!
!!
!!***


#include "Flash.h"
Module sim_interface

  implicit none

  interface
     subroutine sim_ptAnaGetVelComponent(particles,maxParticlesPerProc,&
          varGrid,xPart,yPart,zPart,t,propPart,k)

       integer, INTENT(in) :: maxParticlesPerProc
       real, INTENT(inout),dimension(NPART_PROPS,maxParticlesPerProc) :: particles
       integer, INTENT(in)           :: varGrid, propPart
       integer, INTENT(in)           :: xPart, yPart, zPart, k
       real, INTENT(in)           :: t
     end subroutine sim_ptAnaGetVelComponent
  end interface

  interface
     subroutine sim_getComputedError (error)
       real, INTENT(out)  :: error
     end subroutine sim_getComputedError
  end interface

  interface
     subroutine sim_grAnaGetVelComponent(result,&
          varGrid,xCoord,yCoord,zCoord,t)
       real, INTENT(out)           :: result
       integer, INTENT(in)           :: varGrid
       real, INTENT(in)           :: xCoord, yCoord, zCoord, t
     end subroutine sim_grAnaGetVelComponent
  end interface

end Module sim_interface
