!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/Async/ed_commGetNewRays
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

subroutine ed_commGetNewRays()
  use EnergyDeposition_data, ONLY : ed_rays, ed_rayCount, ed_maxRayCount
  use ed_commData, ONLY : ed_commDebug
  use UTPipeline, ONLY : UTPipeline_getItems
  use ed_commInterface, ONLY : ed_commCheckRayLocation
  implicit none
  integer :: rayCountOld

  rayCountOld = ed_rayCount
  call UTPipeline_getItems(ed_rays, ed_maxRayCount, ed_rayCount)
  if (ed_rayCount > rayCountOld .and. ed_commDebug) then
     call ed_commCheckRayLocation(rayCountOld+1, ed_rayCount)
  end if
end subroutine ed_commGetNewRays
