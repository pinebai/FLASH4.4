!!****if* source/Particles/localAPI/pt_updateTypeDS
!!
!! NAME
!!
!!  pt_updateTypeDS
!!
!! SYNOPSIS
!!
!!  pt_updateTypeDS( integer(in)  :: particlesPerBlk(MAXBLOCKS,NPART_TYPES))
!!
!! DESCRIPTION
!!
!!  Part of the time advancement routine for the Particle Unit.
!!    Calculates the number of each type of particles in the sorted data structure 
!!    (optional return).    Returns the summed total of passive and active particles.
!!    Within the particles data structure, 
!!    PASSIVE_PART_TYPE in located in positions 1:totalPassive
!!    All active particle types are located in positions 
!!             totalPassive+1:totalPassive+1+totalActive
!!
!! ARGUMENTS
!!
!!  particlesPerBlk -- an array of size MAXBLOCKS,NPART_TYPES.  
!!                      It is produced by
!!                      Grid_sortParticles and holds the number 
!!                      of each type of particle in each block.
!!
!! NOTE
!!
!!   The values are actually calculated already in Grid_sortParticles but
!!   encapsulation rules make recalculating there here easier than passing them
!!   back and forth in data, etc.
!!
!!
!!***

!===============================================================================

subroutine pt_updateTypeDS (particlesPerBlk)

  implicit none

#include "Flash.h"
  
  integer, INTENT(in), dimension(MAXBLOCKS,NPART_TYPES) :: particlesPerBlk


  return
  
end subroutine pt_updateTypeDS

!===============================================================================

