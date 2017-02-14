!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commSortParticles
!!
!!  NAME     
!!   ed_commSortParticles
!!
!!  SYNOPSIS
!!   ed_commSortParticles()
!!
!!  DESCRIPTION 
!!    This subroutine sorts the particles in block order.  The underlying
!!    sort happens in Grid_sortParticles or ed_comm_sort_particles.  The
!!    ed_comm_sort_particles subroutine is available if we wish to create
!!    a laser enabled FLASH application without the GridParticles dependency.
!!
!!  ARGUMENTS
!!
!!  NOTES
!!    Rays must be sorted both before and after this subroutine.
!!
!!  SIDE EFFECTS
!!    ed_rays and ed_rayCount will be updated.
!!
!!***

#include "EnergyDeposition.h"

subroutine ed_commSortParticles()
  use EnergyDeposition_data, ONLY : ed_rays, ed_maxRayCount, ed_rayCount
  use ed_commInterface, ONLY : ed_comm_sort_particles
  use Grid_interface, ONLY : Grid_sortParticles
  implicit none
  integer :: raysPerBlk(1:MAXBLOCKS,1:1)
  integer :: initialRayCount
#ifdef USE_GRID_PARTICLES
  logical, parameter :: useGridParticles = .true.
#else
  logical, parameter :: useGridParticles = .false.
#endif

  if (ed_rayCount > 0) then
     if (useGridParticles) then
        call Grid_sortParticles(ed_rays, RAY_ATTR_COUNT, ed_rayCount, 1, &
             ed_maxRayCount, raysPerBlk, RAY_BLCK)
     else
        initialRayCount = ed_rayCount
        call ed_comm_sort_particles(ed_rays, initialRayCount)
        ed_rayCount = count(ed_rays(RAY_BLCK,1:initialRayCount) >= 1.0)
     end if
  end if
end subroutine ed_commSortParticles
