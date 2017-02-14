!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/Simulation_initParticleAttrib
!!
!! NAME
!!
!!  Simulation_initParticleAttrib
!!
!! SYNOPSIS
!!
!!  call Simulation_initParticleAttrib(logical(IN) :: restart)
!!
!! DESCRIPTION
!!
!!  A stub for initializing additional particle attributes if necessary.
!!
!!  This interface is called during initialization after
!!   o  particles positions,
!!   o  particle velocities and other properties needed to advance
!!      particle trajectories,
!!   o  mass properties for active particles,
!!   o  and particle tags
!!  have all been initialized (or possibly read from a checkpoint, if
!!  restart is true).
!!  
!!  This implementation initializes PosInit{X,Y,Z} particle properties
!!  to the corresponding Pos{X,Y,Z} values.
!!
!!
!! ARGUMENTS
!!
!!  restart - true if restarting from a checkpoint, false otherwise.
!!
!! SEE ALSO
!!
!!  Driver_intFlash
!!  Particles_initPositions
!!
!!***

subroutine Simulation_initParticleAttrib(restart)

  use Particles_data, ONLY:  particles, useParticles
  implicit none

#include "Flash.h"

  logical,intent(in) :: restart

  if (.NOT. restart) then
     if (useParticles) then

#ifdef POSINITX_PART_PROP
        particles(POSINITX_PART_PROP,:) = particles(POSX_PART_PROP,:)
#endif
#ifdef POSINITY_PART_PROP
        particles(POSINITY_PART_PROP,:) = particles(POSY_PART_PROP,:)
#endif
#ifdef POSINITZ_PART_PROP
        particles(POSINITZ_PART_PROP,:) = particles(POSZ_PART_PROP,:)
#endif

     end if
  end if

end subroutine Simulation_initParticleAttrib
