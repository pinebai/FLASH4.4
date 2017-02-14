!!****if* source/physics/sourceTerms/EnergyDeposition/EnergyDepositionMain/Laser/LaserComm/ed_commHandleOffBlkRay
!!
!!  NAME     
!!   ed_commHandleOffBlkRay
!!
!!  SYNOPSIS
!!   ed_commHandleOffBlkRay(integer, intent(IN) :: ray)
!!
!!  DESCRIPTION 
!!    This subroutine moves the ray to the correct block.  If the ray exists 
!!    on a block on my MPI rank then it updates the block ID field in ed_rays.
!!    If the ray exists on a block on another MPI rank then the ray is sent
!!    to that MPI rank and ed_rays is set to NONEXISTENT.
!!    
!!  ARGUMENTS
!!    ray : index of the ray in ed_rays
!!
!!  SIDE EFFECTS
!!    ed_rays is updated (in a non-stub implementation).
!!
!!  NOTES
!!    This is a stub version that does not do anything.
!!***

subroutine ed_commHandleOffBlkRay(ray)
  implicit none
  integer, intent(IN) :: ray
end subroutine ed_commHandleOffBlkRay
