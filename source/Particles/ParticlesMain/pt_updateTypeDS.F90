!!****if* source/Particles/ParticlesMain/pt_updateTypeDS
!!
!! NAME
!!
!!  pt_updateTypeDS
!!
!! SYNOPSIS
!!
!!  call pt_updateTypeDS(integer(in) :: particlesPerBlk(MAXBLOCKS,NPART_TYPES))
!!
!! DESCRIPTION
!!
!!  Updates the pt_typeInfo data structure.
!!
!!  Part of the time advancement routine for the Particle Unit.
!!  This routine will usually be called right after a call to Grid_sortParticles.
!!  The particlesPerBlk argument should be the elementsPerBlk returned by
!!  Grid_sortParticles.
!!  Calculates the number of each type of particles in the sorted data structure,
!!  and the index where particles of this type begin in the particles data.
!!
!! EXAMPLE
!!
!!    Assume that there are two classes of particles, 'active' and
!!    passive, in that order.
!!    Then particles of active particle type are located in positions 
!!             totalPassive+1:totalPassive+1+totalActive
!!
!! ARGUMENTS
!!
!!  particlesPerBlk -- an array of size MAXBLOCKS,NPART_TYPES.  It is produced by
!!    Grid_sortParticles and holds the number of each type of particle in each block.
!!
!! SIDE EFFECTS
!!
!!  Updates the pt_typeInfo data structure.
!!
!! NOTES
!!   The values are actually calculated already in Grid_sortParticles but
!!   encapsulation rules make recalculating there here easier than passing them
!!   back and forth in data, etc.
!!
!!   If there is only one particle type, the operation is trivial.
!!   Otherwise, it is assumed that particles are sorted with the particle
!!   type as the sort key of highest priority, as in
!!
!!      call Grid_sortParticles(particles,NPART_PROPS,pt_numLocal,NPART_TYPES, &
!!           pt_maxPerProc,particlesPerBlk,BLK_PART_PROP, TYPE_PART_PROP)
!!
!!
!!***

!===============================================================================

subroutine pt_updateTypeDS (particlesPerBlk)

  use Particles_data, ONLY: pt_typeInfo
  
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Particles.h"
  
  integer, INTENT(in), dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk

  integer :: i,startPoint

  startPoint=1
  do i=1,NPART_TYPES
     pt_typeInfo(PART_TYPE_BEGIN,i)=startPoint
     pt_typeInfo(PART_LOCAL,i) = sum(particlesPerBlk(:,i))
     startPoint=startPoint+pt_typeInfo(PART_LOCAL,i)
  end do

  return
  
end subroutine pt_updateTypeDS

!===============================================================================

