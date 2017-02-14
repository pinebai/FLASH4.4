!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Sync/ed_commGetNewRays
!!
!!  NAME     
!!   ed_commGetNewRays
!!
!!  SYNOPSIS
!!   ed_commGetNewRays()
!!
!!  DESCRIPTION 
!!    This subroutine copies newly received rays from the internal
!!    pipeline buffer into ed_rays.
!!
!!  ARGUMENTS
!!
!!  NOTES
!!    Rays must be sorted both before and after this subroutine.
!!
!!  SIDE EFFECTS
!!    ed_rays and ed_rayCount may be updated.
!!
!!***

#include "EnergyDeposition.h"

subroutine ed_commGetNewRays()
  use Grid_interface, ONLY : Grid_moveParticles
  use EnergyDeposition_data, ONLY : ed_rays, ed_maxRayCount, ed_rayCount, &
       ed_particleIndexList, ed_particleIndexCount, ed_rayDeterminism      
  implicit none

  call Grid_moveParticles (ed_rays,               &
       RAY_ATTR_COUNT,        &
       ed_maxRayCount,        &
       ed_rayCount,           &
       ed_particleIndexList,  &
       ed_particleIndexCount, &
       ed_rayDeterminism      )
end subroutine ed_commGetNewRays
