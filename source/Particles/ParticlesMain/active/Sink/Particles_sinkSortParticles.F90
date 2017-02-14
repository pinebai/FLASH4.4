!!****if* source/Particles/ParticlesMain/active/Sink/Particles_sinkSortParticles
!!
!! NAME
!!
!!  Particles_sinkSortParticles
!!
!! SYNOPSIS
!!
!!  call Particles_sinkSortParticles()
!!
!! DESCRIPTION
!!
!!  Calls Grid_sortParticles().
!!
!! ARGUMENTS
!!
!! NOTES
!!
!!   written by Christoph Federrath, 2014
!!
!!***

subroutine Particles_sinkSortParticles()

  use Particles_sinkData, ONLY : particles_local, localnp,  & 
       sink_maxSinks, pt_sinkParticleProps
  use Grid_interface, ONLY : Grid_sortParticles

  implicit none 

#include "Flash.h"

  integer, dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk
  
  call Grid_sortParticles(particles_local, pt_sinkParticleProps, localnp, &
    NPART_TYPES, sink_maxSinks, particlesPerBlk, BLK_PART_PROP)

end subroutine
