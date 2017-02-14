!!****if* source/Simulation/SimulationMain/unitTest/ParticlesAdvance/HomologousPassive/pt_allBlksInitPos
!!
!! NAME
!!    pt_allBlksInitPos
!!
!! SYNOPSIS
!!    call pt_allBlksInitPos( )
!!
!! DESCRIPTION
!!
!!    Initialize particle positions for all particles.
!!
!!    This default implementation just loops over all local LEAF blocks and does a
!!    separate particle initialization by calling Particles_initPositions for each
!!    of them.
!!
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!    pt_maxPerProc  INTEGER [100]   Maximum number of particles per processor
!!    
!!  NOTES
!!
!!    The particles array should be allocated and initialized to NONEXISTENT before entry.
!!    pt_numLocal should have been set to 0. It will be set to the actual number of
!!    particles initialized on return.
!!
!!  SEE ALSO
!!
!!    Particles_initPositions
!!***

  
subroutine pt_allBlksInitPos ()
  
  use Particles_interface, ONLY : Particles_initPositions
  use Grid_interface, ONLY : Grid_getListOfBlocks
  use Particles_data, ONLY:  particles
  
  implicit none
  
#include "constants.h"
#include "Flash.h"



  integer :: blockCount, blockID
  integer :: blockList(MAXBLOCKS)
  integer :: b
  logical :: posInitialized, updateRefine


!-------------------------------------------------------------------------------

     ! Distribute initial positions
     call Grid_getListOfBlocks(LEAF,blockList,blockCount)

     posInitialized=.false.
     updateRefine=.false.
     call Particles_initPositions(posInitialized, updateRefine)




#ifdef POSINITX_PART_PROP
     particles(POSINITX_PART_PROP,:) = particles(POSX_PART_PROP,:)
#endif
#ifdef POSINITY_PART_PROP
     particles(POSINITY_PART_PROP,:) = particles(POSY_PART_PROP,:)
#endif
#ifdef POSINITZ_PART_PROP
     particles(POSINITZ_PART_PROP,:) = particles(POSZ_PART_PROP,:)
#endif
  
  return

end subroutine pt_allBlksInitPos
